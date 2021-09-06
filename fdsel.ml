module F = Fdsel_lib

(* argument parsing *)
let argspec = ref []
let afusage = ref "Usage:\n  fdsel infer ...\n  fdsel simulate ...\n"
let argsfail msg = (Printf.printf "Failure: %s\n%s" msg
  (Arg.usage_string !argspec !afusage)) ; failwith "Argument failure."
let set_once ref v =
  if !ref = None then ref := Some v
  else raise (Arg.Bad "Only one argument is supported")
let warn str = Printf.eprintf "Warning: %s\n%!" str
let req s z = match !z with Some x -> x | _ -> argsfail s
let def d z = match !z with Some x -> x | _ -> d
let map = Mu.map 
let norm = F.norm

(* parameters from command line are global variables *)
let timeseries_fn = ref None 
and updates_out_fn = ref None 
and cenlev = ref None
and cenub = ref false
and lifespan = ref None 
and params_fn = ref None
and param_mu = ref None
and param_lambda = ref None
and residuals_fn = ref None
and dist_fn = ref None
and rankdist_fn = ref None
and sizes_fn = ref None
and nbins = ref None
and nlogbins = ref None
and nqbins = ref None
and minfreq = ref None
and bootstrap_out_fn = ref None
and bins_in_fn = ref None
and bins_out_fn = ref None
and ne_out_fn = ref None
and run = ref (fun () -> ())
and run_str = ref ""

let breaks_fn = ref None
and ngens = ref None
and pop_size = ref None
and novelty = ref false
and neutrality = ref false
and wmlim = ref false
and demo_fn = ref None
and ts_fn = ref None
and ts_fin = ref false
and dilation = ref None
and contraction = ref None
and burnin = ref None
and seedv = ref None
and nomut = ref false
and nozero = ref false
and wildtype = ref None

let inch fn = if fn = "-" then stdin else open_in fn
let ouch fn = if fn = "-" then stdout else open_out fn

let parse_minfreq pss = match !minfreq, !cenlev with
    (None, None) -> None
  | (Some x, None) -> Some x
  | (None, Some x) -> Some (F.max_censored_freq (x - 1) pss)
  | (Some x, Some y) ->
      warn "Manually overriding censorship minfreq with -f parameter" ;
      warn "This can cause unpredictable results if -f param < cens. minfreq" ;
      Some x

let parse_binning_args minfreq annuds bins_inch =
  match (!nbins, !nlogbins, !nqbins, bins_inch) with
      (None, None, None, None) -> argsfail
        "At least one of -k, -l, -q or -B must be specified."
    | (Some nbins, None, None, None) -> 
        (nbins, F.equal_binning nbins minfreq)
    | (None, Some nbins, None, None) -> 
        (nbins, F.log_binning nbins !wildtype minfreq annuds)
    | (None, None, Some nbins, None) ->
        (nbins, F.quantile_binning nbins !wildtype minfreq annuds)
    | (None, None, None, Some ch) ->
        let breaks = List.sort compare (F.parse_nlsv_ch float_of_string ch) in
        (List.length breaks + 1, (F.explicit_indf breaks, breaks))
    | _ -> argsfail
        "At most one of -k, -l, -q or -B must be specified."

let parse_update_args ts =
  match (!nomut, !nozero, !cenlev, !cenub) with
    (false,false,None,false) -> F.annupdate_data ts
  | (false,false,None,true) -> argsfail "-c argument is required for -C"
  | (false,false,Some _,false) -> F.annupdate_data ts
  | (false,false,Some cenct,true) -> 
      let wtt = req "-w argument is required for -C" wildtype in
      F.annupdates_ub cenct wtt ts
  | (false,true,None,false) -> F.annupdate_data_nozero ts
  | (true,true,None,false) -> warn "-U implies -u" ; F.annupdate_data_nozero ts
  | (true,false,None,false) -> F.annupdate_data_nomut ts
  | _ -> argsfail "-u/-U and -c/-C are incompatible"

let run_inf () =
  let cl = close_out in
  let paramsch = ouch (req "No output file (-o)?" params_fn) in
  let chopt x = Mu.opt_fmap ouch !x in
  let chipt x = Mu.opt_fmap inch !x in
  let doopt oo f = ignore (Mu.opt_fmap f oo) in
  let sizesch = chopt sizes_fn (* open any provided filenames *)
    and bins_inch = chipt bins_in_fn
    and updates_outch = chopt updates_out_fn
    and bootstrap_outch = chopt bootstrap_out_fn 
    and bins_outch = chopt bins_out_fn 
    and residualsch = chopt residuals_fn 
    (* and ne_outch = chopt ne_out_fn *) in
  let ts_gtc = F.parse_timeseries_gtc
	(req "No input timeseries (-i)?" timeseries_fn) in
  let timeseries = match !lifespan with
      None -> F.gen_type_count_to_timeseries ts_gtc
    | Some k -> F.births_gtc_to_timeseries (fun _ -> k) ts_gtc in
  let pop_sizes = F.pop_sizes_ts timeseries in
  let minfreq = parse_minfreq pop_sizes in
  let mutp = match !cenlev with
      Some count -> F.nomut_mutp
    | None -> F.default_mutp in
  let annuds = parse_update_args timeseries in
  doopt sizesch (fun ch -> F.fprint_pop_sizes_ts ch "" timeseries ; cl ch) ;
  let (nbins, (indf, breaks)) = parse_binning_args minfreq annuds bins_inch in
  let (totbins, indty) = F.parse_indty_args nbins indf !wildtype minfreq in
  doopt bins_outch (fun ch -> F.fprint_breaks ch breaks ; cl ch) ;
  let (params, _, dress, lress) =
    F.infer_params !wildtype "MUT" mutp indty totbins annuds in
  doopt updates_outch (fun ch -> F.fprint_diff_ress ch dress; cl ch) ;
  F.fprint_params paramsch params minfreq !wildtype "MUT" mutp indty annuds; 
  cl paramsch ;
  if Mu.opt_some residualsch then (
    let res = F.ts_residuals
      (fst params.F.mu) (F.w_pwc indf (fst params.F.sels)) "MUT" timeseries in
    doopt residualsch (fun ch ->
      F.fprint_ts_residuals ch "" Mu.identity res ; cl ch)) ;
  doopt bootstrap_outch (fun ch ->
    let infer = F.infer_nocis "MUT" mutp indty totbins in
    ignore (F.bootstrap_params_pwc_inc_fprint ch 500 infer annuds); cl ch)
  (* ;
  (if Mu.opt_some residualsch || Mu.opt_some ne_outch then (
    let res = F.ts_residuals
      (fst params.F.mu) (F.w_pwc indf (fst params.F.sels)) "MUT" timeseries in
    doopt residualsch (fun ch ->
      F.fprint_ts_residuals ch "" Mu.identity res ; cl ch) ;
    doopt ne_outch (fun ch ->
      F.fprint_nes ch "" (F.ne_timeseries (F.ts_to_uds_residuals res));cl ch))) *)

let run_sim () =
  let time_dilation = def 1 dilation in
  let pop_contraction = def 1. contraction in
  let contract_ps ps = Mu.round (float_of_int ps *. pop_contraction) in
  let contract_ts pss = map contract_ps pss in
  let params = match !params_fn with
      Some fn -> (match !param_mu with
          Some mu -> { (F.parse_params fn) with F.mu = (mu, 0.) }
          | None -> F.parse_params fn)
    | None -> (if not (!novelty || !neutrality) 
        then argsfail "Parameters file required"
        else {
          F.mu = (req "mu must be specified (-m or -i)" param_mu, 0.);
          F.ne = (1., 1.); F.sels = ([||],[||]) }) in
  let tsoch = ouch (req "Output timeseries is required" timeseries_fn) in
  let breaks = if !novelty || !neutrality then [] 
    else F.parse_nlsv float_of_string 
      (req "Bin breakpoints required" breaks_fn) in
  (if !novelty then (
    (if Array.length (fst params.F.sels) != 1 then 
      warn "ss beyond s1 ignored in novelty simulation") ;
    (if Mu.opt_some !breaks_fn then 
      warn "(-B) Bins ignored in novelty simulation") ;
    (if Mu.opt_some !contraction then
      warn "(-C) contraction unsupported in novelty simulation"))
    else (if !wmlim then warn "-w without -N is meaningless"));
  let tss = Mu.opt_fmap
       (fun fn -> 
         let ts = F.parse_timeseries fn in
         (ts, 
           Mu.sum (map snd (snd (List.hd ts))),
           Mu.sum (map snd (snd (List.hd (Mu.rev ts)))))) 
       !ts_fn in
  let demogv = match (!ngens, !pop_size, !demo_fn, tss) with
      (None, None, None, Some (tsi, _,_)) ->
        F.Variable (contract_ts (F.pop_sizes_ts tsi))
    | (None, None, Some demo, _) ->
        F.Variable (contract_ts (F.parse_nlsv int_of_string demo))
    | (None, None, _, _) ->
        warn "No demography specified; using N=1000 for 100 gen.\n" ;
        F.Fixed (contract_ps 1000, 100)
    | (ng, ps, None, Some (tsi, tsips, tsfps)) ->
        let thisps = Mu.opt_def (if !ts_fin then tsfps else tsips) ps 
        and thisng = Mu.opt_def 100 ng in
        warn (Printf.sprintf "Using fixed demography n=%d g=%d despite -T."
            thisps thisng) ;
        F.Fixed ((contract_ps thisps), thisng)
    | (ng, ps, None, None) ->
        let thisps = Mu.opt_def 1000 ps and thisng = Mu.opt_def 100 ng in
        F.Fixed ((contract_ps thisps), thisng)
    | (Some _, _, Some _, _) -> argsfail "-d and -g are incompatible"
    | (_, Some _, Some _, _) -> argsfail "-d and -n are incompatible" in
  let init_idpop = match (tss, demogv) with
      (Some (ts,_,_), _) -> (* one could contract init pop, but why? *)
        F.idpop_of_pop (snd (List.hd (if !ts_fin then Mu.rev ts else ts)))
    | (None, F.Fixed (ps, _)) -> F.monomorphic_idpop ps
    | (None, F.Variable pss) -> F.monomorphic_idpop (List.hd pss) in
  (* simulate using inferred parameters *)
  let model = if !novelty then (if !wmlim 
      then F.noveltybias_wm_wf (fst params.F.sels).(0)
      else F.noveltybias_wf time_dilation 
        (fst params.F.mu) (fst params.F.sels).(0))
    else (if !neutrality then 
      F.wright_fisher_ time_dilation (Mu.opt_def 1. !param_lambda)
        (fst params.F.mu) F.w_neutral 
    else
      F.wright_fisher_ time_dilation (Mu.opt_def 1. !param_lambda) 
        (fst params.F.mu)
        (F.w_pwc (F.explicit_indf breaks) (fst params.F.sels))) in
  let burnin_init = F.simulate
    ~seed:(Some (Mu.opt_def_f Rand.time_of_day !seedv))
    ~demog:(F.Fixed (Mu.sum (map snd (snd init_idpop)), def 0 burnin))
    F.agg_none (snd init_idpop) model init_idpop in
  (* let burnin_idpop = F.idpop_of_counts 0 (map snd burnin_init) in *)
  let burnin_init = let (mint, maxt) = Mu.min_max (map fst burnin_init) in
    (maxt - mint + 1, map (fun (ty,ct) -> (ty - mint + 1, ct)) burnin_init) in
  F.fprint_pops_tsv_header tsoch "" ;
  F.fprint_pops_tsv_inc string_of_int tsoch "" (0, (snd burnin_init)) ;
  ignore (F.simulate ~seed:(Some (Mu.opt_def_f Rand.time_of_day !seedv))
    ~demog:demogv 
    (F.agg_all_print (F.fprint_pops_tsv_inc string_of_int tsoch "") 1)
    (F.agg_all_knil (snd burnin_init)) model burnin_init)

let run_ts () =
  let cl = close_out in
  let chopt x = Mu.opt_fmap ouch !x in
  let chipt x = Mu.opt_fmap inch !x in
  let optn = Mu.opt_none and opts = Mu.opt_some in
  let doopt oo f = ignore (Mu.opt_fmap f oo) in
  let sizesch = chopt sizes_fn and bins_inch = chipt bins_in_fn
    and updates_outch = chopt updates_out_fn
    and bins_outch = chopt bins_out_fn and distch = chopt dist_fn 
    and rankdistch = chopt rankdist_fn in
  (if optn sizesch && optn bins_outch && optn distch 
    && optn rankdistch && optn updates_outch then failwith "Nothing to do.\n");
  let ts_gtc = F.parse_timeseries_gtc
	(req "No input timeseries (-i)?" timeseries_fn) in
  let timeseries = match !lifespan with
      None -> F.gen_type_count_to_timeseries ts_gtc 
    | Some k -> F.births_gtc_to_timeseries (fun _ -> k) ts_gtc in
  let pop_sizes = F.pop_sizes_ts timeseries in
  let minfreq = parse_minfreq pop_sizes in
  let annuds = parse_update_args timeseries in
  let updates = F.bare_updates annuds in
  doopt updates_outch (fun ch -> F.fprint_updates ch updates) ;
  doopt sizesch (fun ch -> F.fprint_pop_sizes_ts ch "" timeseries ; cl ch) ;
  (if opts bins_outch then (
    let (nbins,(indf,breaks)) = parse_binning_args minfreq annuds bins_inch in
    doopt bins_outch (fun ch -> F.fprint_breaks ch breaks ; cl ch)));
  doopt distch (fun ch -> 
    F.fprint_freq_dist_log ch "" (F.freq_dist_log 50 (map snd timeseries)) ; 
    cl ch) ;
  doopt rankdistch (fun ch -> 
    F.fprint_n_top_rank_freqs ch "" 5000 (map snd timeseries) ; cl ch)

let parse_arguments =
  let usage_v = "Usage:\n  fdsel infer ...\n  fdsel simulate ...\n  fdsel timeseries...\n"  in
  let usage = "Usage: fdsel (infer|simulate|timeseries)" in

  let inf_usage = "Usage:\n  fdsel infer ...\n"  in
  let inf_argspec = [
    ("-i", Arg.String (set_once timeseries_fn), "Input timeseries tsv file");
    ("-P", Arg.String (set_once updates_out_fn), 
      "Output updates datastructure used for inference");
    ("-c", Arg.Int (set_once cenlev), 
      "Censored counts exist below <count> (required for -C)");
    ("-C", Arg.Set cenub, 
      "impute upper bound estimates for censored counts");
    ("-L", Arg.Int (set_once lifespan), 
      "Treat input timeseries as # births that each live -L gens");
    ("-k", Arg.Int (set_once nbins), "(max) Number of evenly-spaced bins to use");
    ("-l", Arg.Int (set_once nlogbins), 
      "(max) Number of log-evenly-spaced bins to use");
    ("-q", Arg.Int (set_once nqbins), "(max) Number of equal-quantile bins");
    ("-f", Arg.Float (set_once minfreq), 
      "Optional minimum frequency for log or quantile binning");
	("-B", Arg.String (set_once bins_in_fn),
      "Read frequency bin boundaries from a list");
    ("-o", Arg.String (set_once params_fn), "Output parameters");
    ("-p", Arg.String (set_once sizes_fn), "Output population sizes");
    ("-S", Arg.String (set_once bootstrap_out_fn), "Output bootstrap");
	("-b", Arg.String (set_once bins_out_fn),
      "Output frequency bin boundaries (breaks)");
    ("-R", Arg.String (set_once residuals_fn), "Output residuals");
    ("-N", Arg.String (set_once ne_out_fn), 
      "Output per-generation N_e estimates");
	("-u", Arg.Set nomut, "Do not impute mutation when a new type appears");
	("-U", Arg.Set nozero, "Do not impute mutation or extinction of any kind");
    ("-w", Arg.String (set_once wildtype), "Consider <type> as wildtype")
  ] in

  let sim_usage = "Usage:\n  fdsel simulate ...\n"  in
  let sim_argspec = [
    ("-i", Arg.String (set_once params_fn), "Parameters (default: neutrality)");
    ("-m", Arg.Float (set_once param_mu), "mu parameter (ignore -i)");
    ("-l", Arg.Float (set_once param_lambda), "lambda parameter (default: 1)");
    ("-B", Arg.String (set_once breaks_fn), 
       "Bin boundaries input file (\\n-separated)");
    ("-o", Arg.String (set_once timeseries_fn), "Output timeseries");
    ("-g", Arg.Int (set_once ngens),
      "Number of generations (default: copy from -d arg, -T arg or 100)");
    ("-n", Arg.Int (set_once pop_size), 
      "Specify fixed population size (default: copy from -T arg or 1000)");
    ("-N", Arg.Set novelty, 
      "Novelty bias model (not exchangeable; ignores -B and -C)");
    ("-K", Arg.Set neutrality, 
      "Neutral model (ignores -B, -C)");
    ("-w", Arg.Set wmlim, "Use weak mutation limit with novelty bias");
    ("-d", Arg.String (set_once demo_fn), "Read population sizes from file");
    ("-T", Arg.String (set_once ts_fn),
      "Copy (initial) population from timeseries (default: monomorphic)");
    ("-F", Arg.Set ts_fin, "(with -T) use final instead of init population");
    ("-D", Arg.Int (set_once dilation), "Time dilation (defult: 1)");
    ("-C", Arg.Float (set_once contraction), 
      "Population contraction (defult: 1.)");
    ("-b", Arg.Int (set_once burnin), "Burn-in generations (default: 0)");
    ("-s", Arg.Int (set_once seedv),
      "Specify RNG seed (default: current time in microseconds)");
  ] in

  let ts_usage = "Usage:\n  fdsel timeseries ...\n"  in
  let ts_argspec = [
    ("-i", Arg.String (set_once timeseries_fn), "Input timeseries tsv file");
    ("-P", Arg.String (set_once updates_out_fn), 
      "Output updates datastructure used for inference");
    ("-c", Arg.Int (set_once cenlev), 
      "Censored counts exist below <count> (required for -C)");
    ("-C", Arg.Set cenub, 
      "impute upper bound estimates for censored counts");
    ("-L", Arg.Int (set_once lifespan), 
      "Treat input timeseries as # births that each live -L gens");
    ("-k", Arg.Int (set_once nbins), "Number of evenly-spaced bins to use");
    ("-l", Arg.Int (set_once nlogbins), 
      "Number of log-evenly-spaced bins to use");
    ("-q", Arg.Int (set_once nqbins), "Number of equal-quantile bins");
    ("-f", Arg.Float (set_once minfreq), 
      "Optional minimum frequency for log or quantile binning");
	("-B", Arg.String (set_once bins_in_fn),
      "Read frequency bin boundaries from a list");
    ("-p", Arg.String (set_once sizes_fn), "Output population sizes");
	("-b", Arg.String (set_once bins_out_fn),
      "Output frequency bin boundaries (breaks)");
    ("-d", Arg.String (set_once dist_fn), "Output frequency distribution");
    ("-Z", Arg.String (set_once rankdist_fn), 
      "Output rank-abundance (Zipf) distribution");
	("-u", Arg.Set nomut, "Do not impute mutation when a new type appears");
	("-U", Arg.Set nozero, "Do not impute any zeros (mutation or extinction)");
  ] in

  let common_usage =
    let noarg arg x = (arg, Arg.Unit (Mu.identity), "") :: x in
    let nohelp x = noarg "-help" (noarg "--help" x) in
    Printf.sprintf 
      "%s\nfdsel infer:%s\nfdsel simulate:%s\nfdsel timeseries:%s\ncommon:" 
      usage_v
      (Arg.usage_string (nohelp inf_argspec) "") 
      (Arg.usage_string (nohelp sim_argspec) "")
      (Arg.usage_string (nohelp ts_argspec) "") in

  let current = ref 0 in
  let anon str = 
    if !current = 1 then match str with
      "simulate" -> run_str := str; argspec := sim_argspec; 
        afusage := sim_usage; run := run_sim
    | "infer" -> run_str := str; argspec := inf_argspec; 
        afusage := inf_usage; run := run_inf
    | "timeseries" -> run_str := str; argspec := ts_argspec; 
        afusage := ts_usage; run := run_ts
    | x -> argsfail (Printf.sprintf "Unmatched argument %s" x) in
  (try Arg.parse_argv_dynamic ~current Sys.argv argspec anon usage
  with Arg.Help message -> Arg.usage [] common_usage ; exit 1
    | Arg.Bad message -> print_string message ; exit 1);
  if not (!run_str = "simulate" || !run_str = "infer" 
    || !run_str = "timeseries") 
  then Arg.usage [] common_usage else ()

let () = !run ()
