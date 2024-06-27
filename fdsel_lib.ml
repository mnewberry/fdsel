let fold = Mu.fold
let map = Mu.map
let iter = Mu.iter
let sumf = Mu.sumf
let cons = Mu.cons
let map2 = Mu.map2
let revh tl h = fold cons tl h
let rev l = revh [] l

let fl = float_of_int
let lnfact n = Gsl.Sf.lnfact n
let norm i z = fl i /. fl z
let normf i z = fl i /. z
let floori x = truncate x
let ceili x = truncate (ceil x)
let atol = Array.to_list
let aofl = Array.of_list
let strcat = String.concat
let strf = string_of_float
let absf = abs_float
let is_nan x = classify_float x = FP_nan

let eq_sels_a xs ys =
  let kons a b acc = acc && (compare a b = 0) in
  Mu.fold2 kons true (atol xs) (atol ys)

let sprint_lof lof =
  Printf.sprintf "[%s]"
    (fold (fun sel st -> Printf.sprintf "%s;%0.3f" st sel)
      (Printf.sprintf "%0.3f" (List.hd lof)) (List.tl lof))

let sprint_mof mof =
  Printf.sprintf "[%s]"
    (fold (fun lof st ->
        let rowst = (fold (fun cell st -> Printf.sprintf "%s;% 18.9f" st cell)
          (Printf.sprintf "% 18.9f" (List.hd (atol lof))) (List.tl (atol lof)))
        in Printf.sprintf "%s%s;\n " st rowst)
      "" (atol mof))

module M = MuMap.Make(struct type t = int let compare = compare end)
module Mp = MuMapP

let sterr n q = 1.96 *. sqrt (q *. (1. -. q) /. fl n)

(* Execute copies of f in nprocs parallel processes and pkons the results
   of each onto pknil.  Once the result of pkons is finishedp, continue
   pkonsing results, but don't spawn any new batches.  f must return an
   integer. *)
let unix_proc_fold nprocs f pkons pknil finishedp =
  let procs = Array.make nprocs
    (None : (int * in_channel * Unix.file_descr * Unix.file_descr) option) in

  let open_proc () =
    let inde, outde = Unix.pipe () in match Unix.fork () with
        0 -> let result = f () in
          (output_binary_int (Unix.out_channel_of_descr outde) result; exit 0)
      | pid -> (pid, Unix.in_channel_of_descr inde, inde, outde) in

  let close_proc _ ch fd1 fd2 =
    let result = input_binary_int ch in
    close_in ch; Unix.close fd2 ; result in

  let rec proc_fold pknil =
    (* If array is filled with none and pknil is finished, we're done *)
    if Array.fold_right (function None -> (&&) true | _ -> (&&) false)
         procs (finishedp pknil) then pknil else
    let new_proc pkn = if finishedp pkn then None else Some (open_proc ()) in
    let kons procn pknil =
      match procs.(procn) with
        None -> (procs.(procn) <- new_proc pknil ; pknil)
      | Some (pid, inch, infd, outfd) -> (
          match Unix.waitpid [Unix.WNOHANG] pid with
              (0, _) -> (ignore (Unix.select [] [] [] 0.05) ; pknil)
            | _ -> let pknil = pkons (close_proc pid inch infd outfd) pknil in
                (procs.(procn) <- new_proc pknil ; pknil)) in
    let pknil = fold kons pknil (Mu.range 0 (nprocs - 1)) in
    proc_fold pknil in
  proc_fold pknil

(* Execute boolean func in batches of batch on nprocs up to max times until
   errt is reached on the fraction of time func returns true.  Returns (total,
   ntrue) of runs.  batch should be high enough that running batch number of
   funcs takes a at least a few seconds to run.  errt is the proportional size
   of the 95% CI. *)
let dogged_monte_carlop label nprocs batch max errt func =
  let finishedp errt (n, ntrue) = if n > max then true else
    let q = fl ntrue /. fl n in
    (n > 0 && ((sterr n q < errt && ntrue < n && sterr n q < errt *. q)
                 || n >= max)) in
  let batchf () = Mu.rec_n batch (fun ct -> if func () then ct+1 else ct) 0 in
  let upf_pkons ntrue (tot, old_ntrue) =
    Printf.printf "UPD-%s\t%d\t%d\n%!" label (tot + batch) (old_ntrue + ntrue);
    (tot + batch, old_ntrue + ntrue) in
  unix_proc_fold nprocs batchf upf_pkons (0, 0) (finishedp errt)

(* {{{ *)

(* Updates datastructure:
   ((from, to) :: rest) :: generations)
   Non-aforeseen types must be registered (0, n). *)

type upd_list = (int * int) list list

(* Note: all three pop_sizes functions differ and give different results.
   However if all types enter and exit the population through zeros (i.e.,
   transitions (0,n) and (n,0) are registered for every type, then the
   following are all equal:
     pop_sizes_ts ts
       = (hd (pop_sizes_uds uds)) :: pop_sizes_uds_prime uds
       = pop_sizes_uds_prime uds @ [hd (rev (pop_sizes_uds_prime uds))]
       = rev (hd (rev (pop_sizes_uds_prime uds)) :: rev (pop_sizes_uds uds))
*)
let pop_size_ts tsgen = fold (fun (_, ct) ps -> ct + ps) 0 tsgen
let pop_sizes_ts timeseries = Mu.map (fun (g,p) -> pop_size_ts p) timeseries

let pop_size ud = fold (fun (ff, _) ss -> ff + ss) 0 ud
let pop_sizes_uds updates = Mu.map pop_size updates

let pop_size_prime ud = fold (fun (_, tt) ss -> tt + ss) 0 ud
let pop_sizes_uds_prime updates = Mu.map pop_size_prime updates

let pop_size_a aud = fold (fun (ty, (ff, _)) ss -> ff + ss) 0 aud
let pop_sizes_auds annups = Mu.map (fun x -> pop_size_a (snd x)) annups

let pop_size_prime_a aud = fold (fun (ty, (_, tt)) ss -> tt + ss) 0 aud
let pop_sizes_auds_prime annups = 
  Mu.map (fun x -> pop_size_prime_a (snd x)) annups

(* Note that print* and fprint* differ in their interpretation of tag *)
let fprint_pop_sizes ouch pre pss =
  Printf.fprintf ouch "%scount\n%!" pre ;
  iter (fun n -> Printf.fprintf ouch "%s%d\n%!" pre n) pss

let print_pop_sizes tag pss = fprint_pop_sizes stdout (tag^"\t") pss

let fprint_pop_sizes_gtc ouch pre gtc =
  let kons (gen,_,cou) mp = Mp.addf gen ((+) cou) 0 mp in
  let mp = fold kons Mp.empty gtc in
  Printf.fprintf ouch "%sgen\tcount\n%!" pre ;
  iter (fun (gen, count) -> 
    Printf.fprintf ouch "%s%d\t%d\n%!" pre gen count) (Mp.bindings mp)

let fprint_pop_sizes_ts ouch pre timeseries =
  Printf.fprintf ouch "%sgen\tcount\n%!" pre ;
  iter (fun (gen, pop) -> 
    let count = pop_size_ts pop in 
    Printf.fprintf ouch "%s%d\t%d\n%!" pre gen count) timeseries

let initial_pop updates =
  fold (fun (frc, toc) (maxtyp, pop) ->
      if frc == 0 then (maxtyp, pop)
      else (maxtyp + 1, (maxtyp + 1, frc) :: pop))
    (0, []) (List.hd updates)

(* update_kons will cons (::) updates onto the updates datastructure (acc)
   unless the type is not found in one or the other population. If either
   population doesn't contain the type, it will either lnfkons (0, ct) or
   tnfkons (ct, 0) onto the update, depending on whether the type was absent in
   this or the previous generation. *)
let update_kons_nf lnfkons tnfkons this_pop (last_pop, acc) =
    let kvm = fold (fun (k, v) m ->
                     BatMap.add k v m) BatMap.empty this_pop in
    let last_kvm = fold (fun (k, v) m ->
                     BatMap.add k v m) BatMap.empty last_pop in
    let kons (ty,last_ct) tl =
      try let this_ct = BatMap.find ty kvm in (last_ct, this_ct) :: tl
      with Not_found -> tnfkons (last_ct, 0) tl in
    let nonmut_ups = fold kons [] last_pop in
    let kons (ty,this_ct) tl =
      try ignore (BatMap.find ty last_kvm); tl
      with Not_found -> lnfkons (0, this_ct) tl in
    (this_pop, (fold kons nonmut_ups this_pop :: acc))

let annupdate_kons_nf lnfkons tnfkons (tgen, tpop) ((lgen, lpop), acc) =
    let tmap = Mp.of_assoc tpop in
    let lmap = Mp.of_assoc lpop in
    let kons (ty,lct) tl =
      try let tct = Mp.find ty tmap in (ty, (lct, tct)) :: tl
      with Not_found -> tnfkons (ty, (lct, 0)) tl in
    let nonmut_ups = fold kons [] lpop in
    let kons (ty,tct) tl =
      try ignore (Mp.find ty lmap); tl
      with Not_found -> lnfkons (ty, (0, tct)) tl in
    ((tgen, tpop), ((lgen, fold kons nonmut_ups tpop) :: acc))

(* impute missing type counts and borrow counts from the wildtype count *)
let annupdate_kons_wt nfct wtt (tgen, tpop) ((lgen, lpop), acc) =
    let tmap = Mp.of_assoc tpop in
    let lmap = Mp.of_assoc lpop in
    let wtc = try Mp.find wtt lmap with Not_found -> 0 in
    let wtpc = try Mp.find wtt tmap with Not_found -> 0 in
    let typecountkons (ty, lct) (wtpc, tl) =
      if ty = wtt then (wtpc, tl) (* don't write update for wt yet *)
      else try let tct = Mp.find ty tmap in 
        (wtpc, (ty, (lct, tct)) :: tl) (* ordinary case *)
      with Not_found -> (* borrow imputed types from wtp count *)
        (wtpc - nfct, (ty, (lct, nfct)) :: tl) in
    let (wtpc, nonmut_ups) = fold typecountkons (wtpc, []) lpop in
    let kons (ty, tct) tl = try ignore (Mp.find ty lmap); tl
      with Not_found -> Mu.cons (ty, (0, tct)) tl in
    ((tgen, tpop), 
      ((lgen, (wtt, (wtc, wtpc)) :: fold kons nonmut_ups tpop) :: acc))

let aggregate_wt wtty indty annuds =
  map (fun (gen, ups) -> 
    let ps = pop_size_a ups in
    let (wtf, wtt, acc) = fold (fun (ty, (frc, toc)) (wtf, wtt, acc) ->
      if indty ps (ty, (frc, toc)) == 0 then (wtf + frc, wtt + toc, acc)
      else (wtf, wtt, (ty, (frc, toc)) :: acc)) (0, 0, []) ups in
    (gen, (wtty, (wtf, wtt)) :: rev acc)) annuds

let annupdates_ub nfct wtt ts =
  rev (snd
    (fold (annupdate_kons_wt nfct wtt) (List.hd ts, []) (List.tl ts)))

let annupdate_kons a b = annupdate_kons_nf Mu.cons Mu.cons a b

let bare_updates annup =
  map (map snd) (map snd annup)

let update_kons a b = update_kons_nf Mu.cons Mu.cons a b

let update_data data =
  rev (snd (fold update_kons ((List.hd data), []) (List.tl data)))
let update_data_k knil data =
  snd (fold update_kons ((List.hd data), knil) (List.tl data))
(* Transitions from zero are not incorporated *)
let update_data_nomut data =
  rev (snd (fold (update_kons_nf (fun a b -> b) Mu.cons)
              ((List.hd data), []) (List.tl data)))
(* Transitions to or from zero are not incorporated *)
let update_data_nozero data =
  rev (snd (fold (update_kons_nf (fun a b -> b) (fun a b -> b))
              ((List.hd data), []) (List.tl data)))

(* See equivalent updates functions above *)
let annupdate_data_nf lnfkons tnfkons timeseries = 
  rev (snd (fold (annupdate_kons_nf lnfkons tnfkons)
              ((List.hd timeseries), []) (List.tl timeseries)))
let idk a b = b (* identity kons means ignore the transition *)
let annupdate_data_nomut ts  = annupdate_data_nf idk     Mu.cons ts
let annupdate_data_nozero ts = annupdate_data_nf idk     idk     ts
let annupdate_data ts        = annupdate_data_nf Mu.cons Mu.cons ts

type 'a pop = ('a * int) list
type idpop = (int * (int pop)) (* maxtype x pop *)

let monomorphic_idpop n = (1, [(1,n)])

let idpop_of_counts maxty cts =
  let (mt, pop) = fold (fun ct (mt, tl) -> if ct = 0 then (mt, tl)
          else (mt + 1, (mt + 1, ct) :: tl)) (maxty, []) cts in
  (mt, rev pop)

let idpop_of_counts_string maxty cts =
  let is = int_of_string and si = string_of_int in
  let (mt, pop) = fold (fun ct (mt, tl) -> if ct = 0 then (mt, tl)
          else (si (is mt + 1), (si (is mt + 1), ct) :: tl)) (maxty, []) cts in
  (mt, rev pop)

(* converts a list of (type : string, count : int) to idpop, attempting to
   convert type strings to integers. If the conversion fails for any string,
   strings are considered opaque and will be numbered sequentially as in
   idpop_of_counts. *)
let idpop_of_pop tcs =
  let kons (ty, ct) (maxty, pop) =
    let tyi = int_of_string ty in
    (max maxty tyi, (tyi, ct) :: pop) in
  try let (mt, pop) = fold kons (0, []) tcs in (mt, rev pop)
  with Failure _ -> idpop_of_counts 0 (map snd tcs)

let pop_of_idpop = snd

let avewv v1 v2 w1 w2 = map2
  (fun x1 x2 -> (w1 *. x1 +. w2 *. x2) /. (w1 +. w2)) v1 v2

let wright_fisher_exp mu wf pop =
  let pimut = mu in
  let ps = pop_size_ts pop in 
  let wis = map (fun (_, ct) -> wf (norm ct ps) *. norm ct ps) pop in
  let sum_wi = sumf wis in
  let piis = map (fun w -> (w /. sum_wi) *. (1. -. pimut)) wis in 
  pimut :: piis

let wright_fisher_exp_ty mu wfty pop =
  let pimut = mu in
  let ps = pop_size_ts pop in 
  let wis = map (fun (ty, ct) -> wfty ps (ty, (ct, 0)) *. norm ct ps) pop in
  let sum_wi = sumf wis in
  let piis = map (fun w -> (w /. sum_wi) *. (1. -. pimut)) wis in 
  pimut :: piis

let neutral_wf ps (maxty, pop) =
  let piis = let ps = pop_size_ts pop in 
    aofl (map (fun (_, ct) -> norm ct ps) pop) in
  let xps = atol (Gsl.Randist.multinomial Rand.rng ~n:ps ~p:piis) in
  (maxty, rev (Mu.fold2 (fun (ty,_) ct pop -> (ty, ct)::pop) [] pop xps))

let wright_fisher_ time_rescaling lambda mu wf ps (maxty, pop) =
  let mult_piis = aofl (wright_fisher_exp mu wf pop) in
  let xps = atol (Gsl.Randist.multinomial Rand.rng ~n:ps ~p:mult_piis) in
  let mut_tot, nm_counts = List.hd xps, List.tl xps in
  let mut_ntypes = if mut_tot > 1 
    then 1 + (Gsl.Randist.binomial Rand.rng ~p:lambda ~n:(mut_tot - 1)) 
    else mut_tot in
  let mut_counts = Array.map ((+) 1) 
    (Gsl.Randist.multinomial Rand.rng ~n:(mut_tot - mut_ntypes)
      ~p:(Array.make mut_ntypes (1. /. fl mut_ntypes))) in
  let maxty, mut_pop = idpop_of_counts maxty (atol mut_counts) in
  let result = (maxty, 
    Mu.fold2 (fun (ty,_) ct pop -> (ty,ct)::pop) mut_pop pop nm_counts) in
  Mu.rec_n (time_rescaling - 1) (neutral_wf ps) result

let wright_fisher lambda mu wf ps (maxty, pop) =
  wright_fisher_ 1 lambda mu wf ps (maxty, pop)

(* The fitness of a type is determined by it's position in the series:
   incrementing type id increments fitness by ss *)
let noveltybias_wf time_rescaling mu ss next_ps (maxty, pop) = 
  let omit_zeros pop = fold 
    (fun (ty, ct) pop -> if ct > 0 then (ty, ct) :: pop else pop) [] pop in
  let this_ps = pop_size_ts pop and pimut = mu in
  let minty = fold (fun (ty, ct) mn -> min ty mn) maxty pop in
  let wis = map (fun (ty, ct) -> exp (fl (ty - minty) *. ss) 
    *. norm ct this_ps) pop in
  let sum_wi = Mu.sumf wis in
  let piis = map (fun w -> (w /. sum_wi) *. (1. -. pimut)) wis in
  let mult_piis = aofl (pimut :: piis) in
  let xps = atol (Gsl.Randist.multinomial Rand.rng ~n:next_ps ~p:mult_piis) in
  let mut_tot, nm_counts = List.hd xps, List.tl xps in
  let maxty, mut_pop = idpop_of_counts maxty (Mu.repeat mut_tot 1) in
  let result = (maxty, omit_zeros
    (Mu.fold2 (fun (ty,_) ct pop -> (ty,ct)::pop) mut_pop pop nm_counts)) in
  Mu.rec_n (time_rescaling - 1) (neutral_wf next_ps) result

(* As above, but simulate "in the weak mutation limit":
   only two types at once are allowed.  Another type is introduced exactly as
   the previous type fixes, skipping the infinite generations of monomorphy.
   Naturally mu is assumed to be infinitesimal. *)
let noveltybias_wm_wf ss next_ps (maxty, pop) = match pop with
    [(ty,_)] -> (maxty + 1, [(maxty+1, 1); (ty, next_ps-1)])
  | _ -> (
    let omit_zeros pop = fold
      (fun (ty, ct) pop -> if ct > 0 then (ty, ct) :: pop else pop) [] pop in
    let this_ps = pop_size_ts pop in
    let minty = fold (fun (ty, ct) mn -> min ty mn) maxty pop in
    let wis = map (fun (ty, ct) -> exp (fl (ty - minty) *. ss) 
      *. norm ct this_ps) pop in
    let sum_wi = Mu.sumf wis in
    let piis = map (fun w -> (w /. sum_wi)) wis in
    let mult_piis = aofl piis in
    let counts = atol (Gsl.Randist.multinomial Rand.rng ~n:next_ps ~p:mult_piis) in
    (maxty, omit_zeros
      (Mu.fold2 (fun (ty,_) ct pop -> (ty,ct)::pop) [] pop counts)))

(* returns (initps, finps, [initfreq], [obsfreq], [expfreq], [neut_expfreq]) *)
let residuals mu wf update =
  let (mutc, ctst, nmctstp1, pst, pstp1) = 
    fold (fun (t1p, t2p) (mc, ctst, ctstp1, pst, pstp1) -> match t1p with
             0->(mc+t2p, ctst, ctstp1, pst, pstp1+t2p)
         | t1p->(mc, t1p::ctst, t2p::ctstp1, pst+t1p, pstp1+t2p))
      (0, [], [], 0, 0) update in
  let (pimut, piis) = 
    let ll = wright_fisher_exp mu wf (snd (idpop_of_counts 0 ctst)) in
    (List.hd ll, List.tl ll) in
  (pst, pstp1,
    (0. :: rev (map (fun ct -> norm ct pst) ctst)),
    (norm mutc pstp1 :: rev (map (fun ct -> norm ct pstp1) nmctstp1)),
    (pimut :: rev piis),
    (pimut :: rev (map (fun ct -> norm ct pst *. (1. -. pimut)) ctst)))

let obsexps mu wf updates =
  Mu.catmap (fun ud -> 
      let (_, _, _, obss, exps, _) = residuals mu wf ud in
      Mu.zip obss exps) updates

let obsexps_of_lhdress ress = 
  fold (fun (gen, psp, types, counts, piis) acc ->
    Mu.fold2 (fun ct pi acc ->
      (norm ct psp, pi) :: acc) acc counts piis) [] ress

let ts_residuals mu wf mutty ts =
  let tcs_to_kv tcs =
    fold (fun (ty,ct) mp -> Mp.add ty ct mp) Mp.empty tcs in
  let add_gen (gen, tcs) ((lastgen, lasttcs), ts_rs) =
    let lookupk (ty, ct) (kv, (mutc, tyups)) =
      (Mp.add ty ct kv,
        try let lastct = Mp.find ty lasttcs in
          (mutc, ((ty, (lastct, ct)) :: tyups))
        with Not_found -> (mutc + ct, tyups)) in
    let (kv, (mutc, tyups)) = fold lookupk (Mp.empty, (0, [])) tcs in
    let (initps, finps, ifs, osfs, expfs, nexpfs) = 
      residuals mu wf ((0, mutc) :: map snd tyups) in
    let nonmut_ty (ty, (fct, tct)) = if fct <> 0 then Some ty else None in
    ((gen, kv),
     (lastgen, initps, finps, 
       (mutty :: Mu.filtmap nonmut_ty tyups), ifs, osfs, expfs, nexpfs) 
       :: ts_rs) in
  let (gen, tcs) = List.hd ts in
  rev (snd (fold add_gen ((gen, tcs_to_kv tcs),[]) (List.tl ts)))

let ts_to_uds_residuals res =
  map (fun (lastgen, initps, finps, tys, ifs, osfs, expfs, nexpfs) ->
    (initps, finps, ifs, osfs, expfs, nexpfs)) res

let ne_estimate obsexps = 
  let ress = aofl (Mu.filtmap (fun (obsf, expf) -> match expf with
            0. -> if expf == 0. then Some 0. else None
          | 1. -> if expf == 1. then Some 0. else None
          | ee -> Some ((obsf -. ee) /. sqrt(ee *. (1. -. ee))))
        obsexps) in
  let var = Gsl.Stats.variance ress in
  let nn = Array.length ress in
  (1. /. var, 
   let x =  2. *. var ** 2. /. fl (nn - 1) in
   x /. var ** 4.)

let ne_timeseries residuals =
  let resss = map (fun (_,_,_,obs,exp,_) -> 
      Mu.filtmap2 (fun obsf expf -> match expf with
            0. -> if expf == 0. then Some 0. else None
          | 1. -> if expf == 1. then Some 0. else None
          | ee -> Some ((obsf -. ee) /. sqrt(ee *. (1. -. ee))))
        obs exp) residuals in
  map (fun ress ->
    let aress = aofl ress in
    let var = Gsl.Stats.variance aress in
    let nn = Array.length aress in 
    (1. /. var, 
     let x =  2. *. var ** 2. /. fl (nn - 1) in
     x /. var ** 4.)) resss

let w_pwc indf sels freq = let ind = indf freq in exp sels.(ind)
let w_pwcty indty sels nn ud = let ind = indty nn ud in exp sels.(ind)
let w_neutral freq = 1.
let wty_neutral _ _ = 1.

let findty indf nnt (ty, (frc, toc)) = indf (norm frc nnt)

let equal_bins_ind minfreq len freq =
  if freq <= 0. || freq > 1.0 then failwith "Frequency outside range"
  else truncate (ceil (freq *. fl len)) - 1

let equal_binning nbins minfreq =
  let mf = Mu.opt_def 1. minfreq in
  (equal_bins_ind mf nbins,
   map (fun i -> mf +. (1. -. mf) *. norm i nbins) (Mu.range 1 (nbins - 1)))

(* Construct some b1 ... bN such that (-inf, b1], (b1, b2], ..., (bN, inf] each
   contain roughly the same number of frequencies in the timeseries.  The b1
   ... bN are unique. If all frequencies in the data are unique, each bin will
   contain at most one more value than any other bin.  However, if there are
   redundant values in the data, the bin contents may differ.
   Ignore frequencies less than minfreq when computing quantiles (hence lower
   bin bounsary is guaranteed to be greater than minfreq). *)
let quantile_breaks wtt nbins minfreq annuds =
  if nbins < 2 then failwith "It's not possible to have less than 2 bins." ;
  let minfreq = Mu.opt_def 0. minfreq in
  let data = (* Convert updates to source frequencies *)
    Mu.catmap (fun (_, uds) -> 
        let tot = pop_size_a uds in
        Mu.catmap (fun (ty, (frc, toc)) -> let ff = norm frc tot in
          if ff >= minfreq && Some ty <> wtt then [ff] else []) uds) annuds in
  let adata = aofl data in
  Array.sort compare adata ; (* sort the source frequencies *)
  let len = Array.length adata in
  (* Don't put in two of the same quantile *)
  let kons ind ((last, skipped), res) =
    let rec next skp n =
      if n = len then failwith (Printf.sprintf 
          "There are too many duplicate frequencies for %n quantile bins." 
          nbins) ;
      if last >= adata.(n)
      then next (skp + 1) (n + 1)
      else ((adata.(n), skp), adata.(n) :: res) in
    next skipped (skipped + ind * (len - 1 - skipped) / nbins) in
  (* skip initial zeros *)
  let nzeros = let rec ct n = if adata.(n) = 0. then ct (n+1) else n in ct 0 in
  let (_, breaks) = fold kons ((0., nzeros), []) 
    (Mu.range 1 (nbins - 1)) in
  rev breaks

let print_data_quantiles tag quantiles =
  Printf.printf "%s\tquantile\n%!" tag ;
  fold (fun x () -> Printf.printf "%s\t%F\n%!" tag x) () quantiles

let fprint_breaks ouch breaks =
  iter (fun x -> Printf.fprintf ouch "%F\n%!" x) breaks

let fprint_nes ouch pre nes =
  Printf.fprintf ouch "%sne\tvar\tlci\tuci\n%!" pre  ;
  iter (fun (ne, ne_var) ->
    let neciw = 1.96 *. sqrt (ne_var) in
    Printf.fprintf ouch "%s%f\t%f\t%f\t%f\n%!" pre 
      ne ne_var (ne-.neciw) (ne+.neciw)) nes

(* Note that there should be one less break than bins because the first first
   and last last bin boundaries are assumed to be at +- inf.  Thus output
   ranges from 0 ... number of breaks, that is, there are nbreaks + 1 = nbins
   possible output values *)
let explicit_indf breaks x =
  let breaks = aofl breaks in
  let seek x n  = if x <= breaks.(n) then n else (n + 1) in
  Mu.rec_n (Array.length breaks) (seek x) 0

let quantile_binning nbins wtt minfreq annuds =
  let breaks = if nbins > 1
  then quantile_breaks wtt nbins minfreq annuds else [] in
  (explicit_indf breaks, breaks)

let log_bins nbins min_freq max_freq =
  let beta = exp ((log max_freq -. log  min_freq) /. fl nbins) in
  let boundaries = if nbins <= 1 then [] else
    Mu.rec_n (nbins - 2) 
      (fun li -> match li with hd :: tl -> hd /. beta :: hd :: tl 
        | [] -> failwith "Can't happen") [max_freq /. beta] in
  (explicit_indf boundaries, boundaries)

let parse_indty_args nbins indf wtty minfreq = match (wtty, minfreq) with
    (Some wtt, Some minf) -> (nbins + 1, 
      (fun nnt (ty, (frc, _)) -> let ff = norm frc nnt in
        if ty = wtt || ff < minf then 0 else indf ff + 1))
  | (Some wtt, None) -> (nbins + 1, 
      (fun nnt (ty, (frc, _)) -> let ff = norm frc nnt in
        if ty = wtt then 0 else indf ff + 1))
  | (None, Some minf) -> (nbins,
      (fun nnt (ty, (frc, _)) -> let ff = norm frc nnt in
        if ff < minf then failwith 
         (Printf.sprintf "No idea how to bin %f = %d/%d with -f." ff frc nnt)
        else indf ff))
  | (None, None) -> (nbins,
      (fun nnt (ty, (frc, _)) -> indf (norm frc nnt)))

let min_max_freq updates =
  let mk cc (mn, mx) ps =
    ((if cc > 0 then min cc mn else mn), max cc mx, ps + cc) in
  fold (fun us (mmn, mmx) ->
      let (mn, mx, ps) = fold (fun (cc,_) (mn, mx, ps) -> mk cc (mn, mx) ps) 
        (max_int, 0, 0) us in
      (min mmn (norm mn ps), max mmx (norm mx ps)))
    (let (mn, mx, ps) = fold (fun (cc,_) (mn, mx, ps) -> mk cc (mn, mx) ps) 
        (max_int, 0, 0) (List.hd updates) in 
      (norm mn ps, norm mx ps))
    (List.tl updates)

let min_max_freq_nowt wtt annupdates =
  fold (fun (_,uds) (mn, mx) ->
    let tot = pop_size_a uds in
    fold (fun (ty, (frc, toc)) (mn, mx) -> let ff = norm frc tot in
      if ff > 0. && wtt <> Some ty then (min ff mn, max ff mx) else (mn, mx))
      (mn, mx) uds) (infinity, neg_infinity) annupdates

let log_binning nbins wtt mfparam updates =
  let (min_freq, max_freq) = min_max_freq_nowt wtt updates in
  log_bins nbins (max (Mu.opt_def 0. mfparam) min_freq) max_freq

(* }}} *)

(* merge_totals pop dist records the marginalized counts
   of the kth most abundant types *)
let merge_totals pop dist =
  let empty = M.empty in
  let insert p map = M.addf p ((+) 1) 0 map in

  let rec merge acc = function
      [], d -> revh d acc
    | (_,p) :: ptl, [] -> merge (insert p empty :: acc) (ptl, [])
    | (_,p) :: ptl, d ::dtl -> merge (insert p d :: acc) (ptl, dtl)
  in merge [] (pop, dist)

(* frequency-based version with n frequency bins *)
let merge_totals_divn n pop dist =
  let tot = Mu.sum (map snd pop) in
  merge_totals (map (fun (a,b) -> (a, (b*2*n + tot)/(2*tot))) pop) dist

let rank_freq_dist_per_n n timeseries =
  fold (merge_totals_divn n) [] timeseries

(* returns [(binid, prob);...] where prob is the proportion of observations
  with frequency in the interval [0,1/nbins),
  [binid/nbins, (binid + 1)/nbins), ..., [(nbins - 1)/nbins, 1] *)

let freq_dist nbins timeseries =
  let insert p map = M.addf p ((+) 1) 0 map in
  let (total, bins) = fold (fun pop (total, bins) ->
      let ps = pop_size_ts pop in
      let bin ct = if ct == ps then nbins - 1 else
        (ct * nbins) / ps in
      fold (fun (_, ct) (total, bins) ->
        (total + 1, insert (bin ct) bins)) (total, bins) pop)
    (0, M.empty) timeseries in
  let (min, _) = M.min_binding bins in
  let (max, _) = M.max_binding bins in
  map (fun id -> (id, norm (try M.find id bins with Not_found -> 0) total))
    (Mu.range min max)

let pop_to_freq timeseries =
  map (fun pop ->
    let ps = pop_size_ts pop in
    map (fun (ty, ct) -> (ty, norm ct ps)) pop) timeseries

let fold_timeseries kons knil timeseries =
  fold (fun pop kn -> fold kons kn pop) knil timeseries

let freq_dist_log nbins timeseries =
  let freqs = pop_to_freq timeseries in
  let (min_freq, max_freq) = let init = snd (List.hd (List.hd freqs)) in
    fold_timeseries (fun y -> Mu.min_max_kons (snd y))
      (init, init) freqs in
  let insert p map = M.addf p ((+) 1) 0 map in
  let (total, bins) = fold (fun pop (total, bins) ->
      let bin freq = if freq == max_freq then nbins - 1 else
        int_of_float ((log freq -. log min_freq) *. fl nbins /.
          (log max_freq -. log min_freq)) in
      fold (fun (_, freq) (total, bins) -> if freq == 0. then (total, bins)
        else (total + 1, insert (bin freq) bins)) (total, bins) pop)
    (0, M.empty) freqs in
  let (min, _) = M.min_binding bins in
  let (max, _) = M.max_binding bins in
  map (fun id ->
        (exp (norm id nbins *. (log max_freq -. log min_freq) +. log min_freq),
        exp (norm (id+1) nbins*.(log max_freq-. log min_freq) +. log min_freq),
        norm (try M.find id bins with Not_found -> 0) total))
    (Mu.range min max)

let clean l = (* TODO: use a sort that behaves well on mostly-sorted lists *)
  (* sort based on population size, largest first *)
  let cmp (a,b) (c,d) = compare (d,c) (b,a) in
  List.sort cmp l 

let idclean (a, l) = (a, List.filter (fun (a, b) -> b <> 0) (clean l))

(* demography is either (popsize, ngens) or a list of popsizes *)
type demography = Fixed of (int * int) | Variable of int list

(* simulates n generations of models of the form (pop -> pop) where a pop is a
   (type * count) list. *)
let simulate ?seed:(seed=None) ?demog:(demog=Fixed (1000, 10)) agg agg_knil 
    model idpop =
  ignore (Mu.opt_fmap Rand.init seed) ;

  let update ps (idpop, aggk) =
    let new_idpop = idclean (model ps idpop) in
    let new_agg = agg (pop_of_idpop new_idpop) aggk in
    (new_idpop, new_agg) in

  match demog with
      Fixed (popsize, ngens) ->
        snd (Mu.rec_n ngens (update popsize) (idpop, agg_knil))
    | Variable pss -> snd (fold update (idpop, agg_knil) pss)

let agg_none x y = x
let agg_none_knil = []

type 'a agg_t = {
  pop : int * (('a * int) list) ;
  ups : (int * (int * (int * int)) list) list ;
  dist : (int M.t) list ;
  pops : (int * ('a * int) list) list ; (* [(generation, [(type,count)])] *)
  nupd : int }



let agg_all_print print sampling_interval pop agg =
  if (agg.nupd + 1) mod sampling_interval <> 0 then {
    agg with nupd = agg.nupd + 1 } else
  let this_pop = agg.nupd + 1, pop in
  print this_pop ;
  { pop = this_pop ;
    ups = snd (annupdate_kons this_pop (agg.pop, agg.ups)) ;
    dist = merge_totals pop agg.dist ;
    pops = this_pop :: agg.pops ;
    nupd = agg.nupd + 1 }

let agg_all sampling_interval pop agg =
  agg_all_print (fun _ -> ()) sampling_interval pop agg

let agg_all_knil pop = {pop=(0,pop); ups=[]; dist=[]; pops=[(0,pop)]; nupd=0}

let fold_updates kons knil updates =
  fold (fun ups kn -> fold kons kn ups) knil updates

let fold_annupdates kons knil updates =
  fold (fun (gen, ups) kn -> fold (kons gen) kn ups) knil updates

(* Pretty-printers {{{ *)

(* 'dist' throughout refers to the rank-abundance distribution. *)
let print_dist n dist =
  let rec pd listpos = function
      [] -> ()
    | hd :: tl ->
      (* Printf.printf "UPD: %d:\nUPD:" listpos ; *)
      fold (fun (k, v) () -> Printf.printf " %d:%0.5f" k (norm v n )) ()
        (M.bindings hd);
      Printf.printf "\n%!" ;
      pd (listpos + 1) tl
   in
   Printf.printf "\x1B\x63" ;
   pd 1 dist

let print_dist_tsv tag n dist =
  let rec pd lpos = function
      [] -> ()
    | hd :: tl ->
      fold (fun (k, v) () ->
        Printf.printf "%s\t%d\t%d\t%0.5f\n" tag lpos k (norm v n))
        () (M.bindings hd);
      pd (lpos + 1) tl
   in
   Printf.printf "%s\trank\tabund\tprob\n" tag ;
   pd 1 dist

let fprint_n_top_rank_freqs ouch pre n timeseries =
  Printf.fprintf ouch "%srank\tfreq\n" pre ;
  iter
    (fun pop ->
      let ps = pop_size_ts pop in
      ignore (Mu.rec_n n (fun (rank, rest) -> match rest with
          [] -> (0, [])
        | (tp,ct) :: tl ->
            Printf.fprintf ouch "%s%d\t%F\n" pre rank (norm ct ps) ;
            (rank + 1, tl)) (1, pop)))
    timeseries

let print_n_top_rank_freqs tag n timeseries =
  fprint_n_top_rank_freqs stdout (Printf.sprintf "%s\t" tag) n timeseries

let print_freq_dist tag nbins dist =
  Printf.printf "%s\tnbins\tbin\tprob\n" tag ;
  iter (fun (bin, freq) ->
    Printf.printf "%s\t%d\t%d\t%F\n" tag nbins bin freq) dist

let fprint_freq_dist_log ch pre dist =
  Printf.fprintf ch "%sbmin\tbmax\tprob\n" pre ;
  iter (fun (min, max, freq) ->
    Printf.fprintf ch "%s%F\t%F\t%F\n" pre min max freq) dist

let print_freq_dist_log tag dist =
  fprint_freq_dist_log stdout (Printf.sprintf "%s\t" tag) dist

let print_pop pop =
  Printf.printf "[" ;
  fold (fun (n,p) () -> Printf.printf "%dx%s;" p n) () pop ;
  Printf.printf "[]]\n"

let print_pop_tsv string_of_type prefix pop =
  Printf.printf "%s\ttype\tcount\n%!" prefix ;
  iter (fun (ty,ct) -> Printf.printf "%s\t%d\t%d\n%!" prefix ty ct) pop

let fprint_pops_tsv string_of_type ouch pre pops =
  Printf.fprintf ouch "%sgen\ttype\tcount\n%!" pre ;
  iter (fun (gen, pop) ->
    iter (fun (ty,ct) ->
      Printf.fprintf ouch "%s%d\t%s\t%d\n%!"
        pre gen (string_of_type ty) ct) pop) pops

let fprint_pops_tsv_header ouch pre =
  Printf.fprintf ouch "%sgen\ttype\tcount\n%!" pre

let fprint_pops_tsv_inc string_of_type ouch pre (gen, pop) =
  iter (fun (ty,ct) ->
    Printf.fprintf ouch "%s%d\t%s\t%d\n%!"
      pre gen (string_of_type ty) ct) pop 

let print_pops_tsv string_of_type tag pops =
  fprint_pops_tsv string_of_type stdout (tag^"\t") pops 

let print_timeseries = print_pops_tsv

let print_updates_dbg tag ups =
  iter (fun ups -> Printf.printf "updates_dbg\t%s\t%s\n%!" tag
    (fold (fun (l,t) acc -> Printf.sprintf "%s%d->%d;" acc l t) "" ups))
    ups

(* unlike the above, print in an ocaml-readable datastructure *)
let fprint_updates ouch updates =
  Printf.fprintf ouch "let updates = [\n%!" ;
  iter (fun ups ->
    Printf.fprintf ouch "  [%s];\n%!"
      (String.concat ";" (map (fun (l,t) ->
        Printf.sprintf "(%d,%d)" l t) ups)))
    updates ;
  Printf.fprintf ouch "] ;;\n"

let print_wp tag ps wp =
  Printf.printf "%s\tct\twp\n%!" tag ;
  let pr ct = Printf.printf "%s\t%d\t%F\n%!" tag ct (wp ct) ; ct + 1 in
  ignore (Mu.rec_n (ps + 1) pr 0)

let fprint_bindata ouch prefix indf breaks updates =
  let breaks = Array.of_list breaks in
  let dd = Array.length breaks + 1 in
  let (bin_gens, bin_transs, bin_fsum, bin_lnfsum) =
    let bin_gens = Array.make dd 0 in
    let bin_transs = Array.make dd 0 in
    let bin_fsum = Array.make dd 0. in
    let bin_lnfsum = Array.make dd 0. in
    iter
      (fun ups ->
        let bin_counts_gen = Array.make dd 0 in
        let nnt = pop_size ups in
        let ff (frc, toc) =
          if frc <> 0 then 
            let freq = norm frc nnt in
            let ind = indf freq in (* ind = b(x_j^t/N_t) *)
            bin_counts_gen.(ind) <- 1 ;
            bin_transs.(ind) <- bin_transs.(ind) + 1 ;
            bin_fsum.(ind) <- bin_fsum.(ind) +. freq ;
            bin_lnfsum.(ind) <- bin_lnfsum.(ind) +. log freq in
        iter ff ups ;
        Array.iteri
          (fun ind bin_ct -> bin_gens.(ind) <- bin_gens.(ind) + bin_ct)
          bin_counts_gen)
      updates ;
    (bin_gens, bin_transs, bin_fsum, bin_lnfsum) in
  let nupd = Mu.sum (map List.length updates) in
  let nnmupd = Mu.sum (map List.length 
    (map (List.filter (fun (x,_) -> x > 0)) updates)) in
  let bounds ind = 
    ((if ind - 1 < 0 then neg_infinity else breaks.(ind - 1)),
     (if ind > dd - 2 then infinity else breaks.(ind))) in

  Printf.fprintf ouch 
    "%snupd\tnnmuupd\tngens\tbin\tllim\tulim\tgenc\ttransc\tmeanfreq\tgeomeanfreq\n%!" prefix ; 
  Mu.iter (fun ind ->
      Printf.fprintf ouch "%s%d\t%d\t%d\t%d\t%e\t%e\t%d\t%d\t%F\t%F\n%!" prefix
        nupd nnmupd (List.length updates)
        ind (fst (bounds ind)) (snd (bounds ind))
        bin_gens.(ind) bin_transs.(ind)
        (bin_fsum.(ind) /. float_of_int bin_transs.(ind))
        (exp (bin_lnfsum.(ind) /. float_of_int bin_transs.(ind))))
    (Mu.range 0 (dd - 1))

(* }}} *)

(* MLEst {{{ *)
let likelihood mu wf updates =
  let sum_gen us tot =
    let (mut_us, nonmut_us) = Mu.partition (fun (a,_) -> a == 0) us in
    let oldps = Mu.sum (map fst us) in
    let mut_ct = Mu.sum (map snd mut_us) in
    (* let nonmut_ct = Mu.sum (map snd nonmut_us) in *)
    (* let newps = mut_ct + nonmut_ct in *)
    let pimut = mu in
    let wis = map (fun (frc,_) -> wf (norm frc oldps) *. norm frc oldps)
                nonmut_us in
    let tot_wi = Mu.sumf wis in
    let piis = pimut :: (map (fun wi -> wi /. tot_wi *. (1. -. pimut)) wis) in
    let counts = mut_ct :: map snd nonmut_us in
    (* Printf.printf "gslll\t%F\t%F\n%!" (Mu.sumf piis) tot_wi ; *)
    tot +. Gsl.Randist.multinomial_lnpdf ~p:(aofl piis) ~n:(aofl counts) in
  fold sum_gen 0. updates

let ml_mutation_rate updates = 
  let (num, denom) = fold_updates
    (fun (frc, toc) (nsum, dsum) ->
      if frc == 0 then (nsum + toc, dsum + toc) else (nsum, dsum + toc))
    (0, 0) updates in
  norm num denom

(* Derivative of the piecewise constant loglhd function wrt all dd parameters.
  The derivative of the log likelihood can be thought of as 'obs = expect'.
  The return value is two arrays, where obs.(k) and expect.(k) are those
  components of the derivative with respect to k. The arguments coef and sels
  are both vectors of parameters. The difference is that in mm, coeff are sk
  that will become the coefficients e^sk to the RHS numerator of Eq 13, that
  is, the variable parameters with respect to the optimization, and sels are
  the iterates or the sk that appear in the RHS denominator. The derivative
  with respect to sk at a given set of parameters sels is:
    let (obs, ex) = diffll sels sels updates in obs.(k) - ex.(k)
  The optimum s_k in the context of estimated sels sms is:
    let (obs, ex) = diffll (Array.map (fun _ -> 0.) sms) sms updates in
      log (obss.(k) /. ex.(k))
  This is possible because exp 0. = 1., and coef is named that to suggest that
  those variable sks only show up as the e^sk coefficient to the ex.  Passing
  sk = 0 in coeff allows you to solve for e^sk as (obs.(k)/ex.(k)) by replacing
  the factor of e^sk with 1. *)
let diffll coef indf sels updates =
    let dd = Array.length sels in (* number of bins *)
    let obs = Array.make dd 0. in (* accumulator for Eq 13 LHS by bin *)
    let exs = Array.make dd 0. in (* acc for Eq 13 RHS by bin *)
    iter (fun ups ->
      let bin_sums = Array.make dd 0. in
      let nnt = pop_size ups in
      let (nonmutc, zz) =
        let kons (frc, toc) (nonmutc, zz) =
          let frcexpsk =
            if frc <> 0 then
              (let ind = indf (norm frc nnt) in
                obs.(ind) <- obs.(ind) +. fl toc;
                bin_sums.(ind) <- bin_sums.(ind) +. fl frc ;
                exp sels.(ind) *. (fl frc))
            else 0. in
          ((if frc = 0 then nonmutc else nonmutc + toc), zz +. frcexpsk) in
        fold kons (0, 0.) ups in
      Array.iteri
        (fun ind bin_sum ->
          exs.(ind) <- exs.(ind) +.
            exp coef.(ind) *. bin_sum *. fl nonmutc /. zz)
        bin_sums) updates ;
    (obs, exs)

(* similar to above, but takes annotated update and mutp, a predicate
   (typically (frc == 0)) to determine whether that transition is a mutation *)
let default_mutp _ (ty, (frc, toc)) = frc = 0
let nomut_mutp _ _ = false
let count_mutp count _ (ty, (frc, toc)) = frc < count
let freq_mutp freq nnt (ty, (frc, toc)) = norm frc nnt < freq

(* mutwt -- use mutp to determine wildtype in addition to wt. *)
let diffll_ress coef sels indty mutp annupds =
  let dd = Array.length sels in (* number of bins *)
  let obs = Array.make dd 0 in (* accumulator for Eq 13 LHS by bin *)
  let exs = Array.make dd 0. in (* acc for Eq 13 RHS by bin *)
  let ress = map (fun (gen, ups) ->
    let bin_sums = Array.make dd 0. in
    let nnt = pop_size_a ups in
    let (mutc, nonmutc, zz, tyress) =
      let kons ((ty, (frc, toc)) as trans) (mutc, nonmutc, zz, tyracc) =
        let ind = indty nnt (ty, (frc, toc)) in
        let frcexpsk =
          if mutp nnt trans then 0.
          else (obs.(ind) <- obs.(ind) + toc;
             bin_sums.(ind) <- bin_sums.(ind) +. fl frc ;
             exp sels.(ind) *. (fl frc)) in
        ((if mutp nnt trans then mutc + toc else mutc),
         (if mutp nnt trans then nonmutc else nonmutc + toc), 
         zz +. frcexpsk, 
         (ty, ind, frc, toc, frcexpsk) :: tyracc) in
      fold kons (0, 0, 0., []) ups in
    Array.iteri
      (fun ind bin_sum ->
        exs.(ind) <- exs.(ind) +.
          exp coef.(ind) *. bin_sum *. fl nonmutc /. zz)
      bin_sums ;
    (gen, nnt, mutc, nonmutc, zz, tyress)) annupds in
  (obs, exs, ress)

(* Legacy version, does not deal with censorship *)
let ml_mm_pwclf logchopt indf sels updates =
  let dd = Array.length sels in

  (* === logging mm iterates to a file, if desired === *)
  ignore (Mu.opt_fmap
    (fun logch ->  Printf.fprintf logch "ml_mm_pwclf\tll\t%s\n%!"
      (strcat "\t" (map (Printf.sprintf "s%d") (Mu.range 1 dd)))) 
    logchopt) ;
  let last_likelihood = ref (-1./.0.) in
  let print_sels lhd sels =
    ignore (Mu.opt_fmap 
      (fun logch ->
        Printf.fprintf logch "ml_mm_pwclf\t%f\t%s\n%!"
          lhd (strcat "\t" (map strf sels)))
      logchopt) ;
    last_likelihood := lhd in
  (* === end logging === *)

  let muest = ml_mutation_rate updates in

  let mm sms =
    let (obs, exs) = diffll (Array.make dd 0.) indf sms updates in
    let optss = Array.mapi (fun ind _ -> log (obs.(ind) /. exs.(ind))) obs in
    let smean = Mu.meanf_n (atol optss) in
    let optss = Array.map (fun x -> x -. smean) optss in
    optss in

  (* recurse mm to fixed likelihood *)
  let rec fix fn init init_lhd =
    let next = fn init in
    let lhd = likelihood muest (w_pwc indf next) updates in
    print_sels lhd (atol next) ;
    if eq_sels_a init next || lhd <= init_lhd then (init, init_lhd)
    else fix fn next lhd in

  let (sels, ll) = fix mm sels neg_infinity in
  (muest, sels, ll)

let annupdates_residuals sels indty mutp annups =
  let dd = Array.length sels in let zeros = Array.make dd 0. in
  let (_, _, ress) = diffll_ress zeros sels indty mutp annups in
  ress

let fprint_diff_ress ouch ress =
  Printf.fprintf ouch
"gen\tinit_ps\tmutc\tnonmutc\tfin_ps\tdenom\ttype\tind\tinitc\tfinc\tinitf\tobsf\tcesk\tnexpf\texpf\trs_inc\n" ;
  iter (fun (gen, nnt, mutc, nonmutc, zz, tyress) ->
    iter (fun (ty, ind, frc, toc, frcexpsk) ->
            let initf = norm frc nnt and obsf = norm toc (mutc + nonmutc) in
            let nonmutf = norm nonmutc (mutc + nonmutc) in
            let nexpf = initf *. nonmutf in
            let expf = frcexpsk /. zz *. nonmutf in
            let rs_inc = (obsf -. expf) /. sqrt (expf *. (1. -. expf)) in
            Printf.fprintf ouch 
"%d\t%d\t%d\t%d\t%d\t%F\t%s\t%d\t%d\t%d\t%F\t%F\t%F\t%F\t%F\t%F\n"
              gen nnt mutc nonmutc (mutc + nonmutc)
              zz ty ind frc toc 
              initf obsf
              frcexpsk
              nexpf expf rs_inc) 
         tyress) 
    ress

(* confidence intervals *)

(* A brief intro on gsl's matrix operations, we write a transformation matrix
   to rotate a vector by 60 degrees around the origin, then we invert the
   matrix and rotate it back:
  let pi = acos (-1.) ;;
  let pio3 = pi /. 3. ;;
  let rot = Gsl.Matrix.of_arrays [|[| cos pio3 ; -. sin pio3 |];
                                   [| sin pio3 ;    cos pio3 |]|] ;;
  rot.{0,0} ;;
  rot.{1,0} ;;
  rot.{0,1} ;;
  let xhat = Gsl.Vector.of_array [|1.;0.|] ;;
  let vzero = let v = Gsl.Vector.copy xhat in Gsl.Vector.set_zero v ; v ;;
  Gsl.Vector.to_array vzero ;;
  Gsl.Blas.gemv Gsl.Blas.NoTrans 1. rot xhat 0. vzero ;;
  Gsl.Vector.to_array vzero ;;
  let rotinv = match Gsl.Linalg.invert_LU (`M rot)
    with `M x -> x | _ -> raise Not_found ;;
  let vzero2 = let v = Gsl.Vector.copy xhat in Gsl.Vector.set_zero v ; v ;;
  Gsl.Blas.gemv Gsl.Blas.NoTrans 1. rotinv vzero 0. vzero2 ;;
  Gsl.Vector.to_array vzero2 ;; *)
let observed_information_var mlmu indty mlss annupdates =
  let dd = Array.length mlss in

  (* complicated index manipulation to divert matrix around nans *)
  let (ddd, ind, inds, indsinv) = 
    Array.fold_left (fun (ddd, ind, inds, indsinv) ss ->
      if is_nan ss then (ddd, ind + 1, (-1) :: inds, indsinv) 
      else (ddd+1, ind+1, ddd :: inds, ind :: indsinv))
    (0, 0, [], []) mlss in
  let inds = aofl (rev inds) in
  let indsinv = aofl (rev indsinv) in

  (* print_endline (sprint_lof (map fl (atol inds))) ;
  print_endline (sprint_lof (atol mlss)) ;
  print_endline (sprint_lof (map fl (atol indsinv))) ; *)

  let varmu =
    let totps = fl (Mu.sum (pop_sizes_auds_prime annupdates)) in
    let (nons, muts) = ((1. -. mlmu) *. totps, mlmu *. totps) in
    let obs_inf = nons /. ((1. -. mlmu)**2.) +. muts /. (mlmu**2.) in
    1. /. obs_inf in

  if ddd < 2 then (varmu, Array.make dd nan) else (

  let bin_sums_by_gen =
    let sums (_, ups) =
      let bin_sums = Array.make ddd 0. in
      let nnt = pop_size_a ups in
      let iter_kons (ty, (frc, toc)) nonmutc =
        if frc = 0 then nonmutc
        else (let ind = indty nnt (ty, (frc, toc)) in
          (* Printf.printf "setting bin_sums.(%d) to %F\n%!" ind
               (bin_sums.(ind) +. exp mlss.(ind) *. fl frc); *)
          (* Printf.printf "%d->%d\n" ind (inds.(ind)) ; *)
          bin_sums.(inds.(ind)) <- bin_sums.(inds.(ind))
            +. exp mlss.(ind) *. fl frc ;
          nonmutc + toc) in
      let nonmutc = fold iter_kons 0 ups in
      let (a1, a2) = (nonmutc, atol bin_sums) in
      (* Printf.printf "gen nonmutc=%d sums=%s\n%!" a1 (sprint_lof a2) ; *)
      (a1,a2)
    in map sums annupdates in

  let d2ll ia ib =
    let kons (nonmutc, bin_sums) tot =
      let (ks, ls, ts) = Mu.fold2 (fun vv ind (ks, ls, ts) ->
        ((if ind = ia then ks +. vv else ks),
         (if ia = ib then (if ind = ia then ls else ls +. vv)
            else (if ind = ib then ls +. vv else ls)),
         (ts +. vv))) (0., 0., 0.) bin_sums (Mu.range 0 (ddd - 1)) in
      tot +. fl nonmutc *. ks *. ls /. (ts ** 2.) *.
               (if ia = ib then 1. else -1.) in
    fold kons 0. bin_sums_by_gen in

  let lower_partials = Gsl.Matrix.create ~init:0. (ddd - 1) (ddd - 1) in
  (* Compute the upper triangle of the observed information matrix,
     set lower to upper *)
  let di = ddd - 1 in
  for i = 0 to di - 1 do for j = 0 to i do
    if i = j then
      lower_partials.{i,j} <- d2ll i i -. 2. *. d2ll i di +. d2ll di di
    else begin
      let d2llpdidj = d2ll i j -. d2ll i di -. d2ll j di +. d2ll di di in
      lower_partials.{i,j} <- d2llpdidj ;
      lower_partials.{j,i} <- d2llpdidj
  end done done ;

  let lower_var_ss = match Gsl.Linalg.invert_LU (`M lower_partials)
    with `M x -> x | _ -> failwith "Could not invert partials" in

  let sdvect = Gsl.Vector.of_array (Array.make (ddd - 1) (-.1.)) in
  let varsd =
    let v1 = Gsl.Vector.copy sdvect in
    Gsl.Blas.gemv Gsl.Blas.NoTrans 
      ~alpha:1. ~a:lower_var_ss ~x:sdvect ~beta:0.  ~y:v1 ;
    Gsl.Blas.dot v1 sdvect in

  let varss = Array.make dd nan in
  varss.(indsinv.(ddd - 1)) <- varsd ;
  iter (fun i -> varss.(indsinv.(i)) <- lower_var_ss.{i,i}) 
    (Mu.range 0 (ddd - 2)) ;
  (varmu, varss))

  (* (* check that setting s_1 to \sum_{i = 2}^D gives the same result as
        setting s_D to \sum_{i = 1}^{D - 1} *)
  let upper_partials = Gsl.Matrix.create ~init:0. (dd - 1) (dd - 1) in
  (* Compute the upper triangle of the observed information matrix,
     set lower to upper *)
  let di = 0 in
  for i = 1 to dd - 1 do for j = 1 to i do
    if i = j then
      upper_partials.{i-1,j-1} <- d2ll i i -. 2. *. d2ll i di +. d2ll di di
    else begin
      let d2llpdidj = d2ll i j -. d2ll i di -. d2ll j di +. d2ll di di in
      upper_partials.{i-1,j-1} <- d2llpdidj ;
      upper_partials.{j-1,i-1} <- d2llpdidj
  end done done ;

  Printf.printf "Observed Information Matrix for s_2...s_D:\n%s\n%!"
    (sprint_mof (Gsl.Matrix.to_arrays upper_partials)) ;
  let upper_det = Gsl.Linalg.det_LU (`M upper_partials) in
  Printf.printf "determinant\t%F\n%!" upper_det ;
  let upper_var_ss = match Gsl.Linalg.invert_LU (`M upper_partials)
    with `M x -> x | _ -> failwith "Could not invert partials" in
  Printf.printf "Inverse of OI Matrix:\n%s\n%!"
    (sprint_mof (Gsl.Matrix.to_arrays upper_var_ss)) ;

  let sdvect = Gsl.Vector.of_array (Array.make (dd - 1) (-.1.)) in
  let varsd =
    let v1 = Gsl.Vector.copy sdvect in
    Gsl.Blas.gemv Gsl.Blas.NoTrans 1. upper_var_ss sdvect 0. v1 ;
    Gsl.Blas.dot v1 sdvect in
  Printf.printf "var(s_1) = %.9f\n%!" varsd *)

  (* (* check partial derivatives by finite difference *)
  let f e1 e2 =
    let sels = Array.copy mlss in
    sels.(0) <- sels.(0) +. e1 ;
    sels.(1) <- sels.(1) +. e2 ;
    likelihood mlmu (w_pwc indf sels) updates in
  let eps = 0.001 and zz = 0. in
  let d00 = (f eps zz -. 2. *. f zz zz +. f (-.eps) zz) /. (eps ** 2.) in
  let d11 = (f zz eps -. 2. *. f zz zz +. f zz (-.eps)) /. (eps ** 2.) in
  let d01 = (f eps eps -. f eps (-.eps) -. f (-.eps) eps +. f (-.eps) (-.eps))
    /. (4. *. eps ** 2.) in

  Printf.printf "d00 %F\n%!" d00 ;
  Printf.printf "d11 %F\n%!" d11 ;
  Printf.printf "d01 %F\n%!" d01 ; *)

let vars_to_cis ms vs =
  let dd = Array.length ms in
  let (l,u) = (Array.make dd 0., Array.make dd 0.) in
  Array.iteri (fun i v -> let w = 1.96 *. sqrt v in
    l.(i) <- ms.(i) -. w ;
    u.(i) <- ms.(i) +. w) vs ; (l, u)

let replacement_fitness mu indf sels updates =
  let kons ups (ftot, ngen) =
      let nnt = pop_size ups in
      let wsum = sumf
        (map (fun (x,_) ->
          if x = 0 then 0.
          else fl x *. exp sels.(indf (norm x nnt))) ups) in
      (ftot +. log (wsum /. fl nnt +. mu), ngen + 1) in
  let (ftot, ngen) = fold kons (0., 0) updates in
  ftot /. fl ngen

(* resurrected code from the past; be skeptical *)
let mean_fitness indf sels updates =
  let dd = Array.length sels in
  let bin_counts = Array.make dd 0 in
  let tot = fold (fun ups n ->
      let sum = Mu.sum (map fst ups) in
      fold (fun (frc, _) n -> if frc = 0 then n else
          (let ind = indf (norm frc sum) in
          bin_counts.(ind) <- bin_counts.(ind) + frc ;
          n + frc))
        n ups)
     0 updates in
  log (fold (fun i sum -> sum +. norm bin_counts.(i) tot *. exp sels.(i))
         0. (Mu.range 0 (dd - 1)))

(* matters of censorship and imputation *)

(* cencount + 0.5 / (min population size) *)
let max_censored_freq cencount pss = 
  Mu.max (map (norm (cencount * 2 + 1)) (map (( * ) 2) pss))

let ml_numigration_rate mutp annupds =
  let (nsum, dsum) =
    fold (fun (_,gen) acc ->
        let nnt = pop_size_a gen in
        fold (fun ((ty,(frc,toc)) as trans) (nsum, dsum) ->
            if mutp nnt trans
            then (nsum + toc, dsum + toc) 
            else (nsum, dsum + toc))
          acc gen)
      (0, 0) annupds in
  norm nsum dsum

let num_replacement_fitness nu indty sels mutp annuds =
  let kons (gen, aups) (ftot, ngen) =
      let nnt = pop_size_a aups in
      let wsum = sumf
        (map (fun ((ty, (x,_)) as trans) ->
          if mutp nnt trans then 0.
          else fl x *. exp sels.(indty nnt trans)) aups) in
      (ftot +. log (wsum /. fl nnt +. nu), ngen + 1) in
  let (ftot, ngen) = fold kons (0., 0) annuds in
  ftot /. fl ngen

(* returns (ll, ress) *)
let likelihood_cen nu wty mutty mutp annupdates =
  let sum_gen (gen, us) (tot, ress) =
    let oldps = pop_size_a us in
    let (mut_us, nonmut_us) = Mu.partition (mutp oldps) us in
    let mut_ct = pop_size_prime_a mut_us in
    (* let nonmut_ct = Mu.sum (map snd nonmut_us) in *)
    (* let newps = mut_ct + nonmut_ct in *)
    let pi0 = nu in
    let (types, counts, nmps, wis, zz) =
      fold (fun (ty,(frc,toc)) (types, cts, nmps, wis, zz) ->
        let wi = wty oldps (ty,(frc,toc)) *. fl frc in
        (ty :: types, toc :: cts, toc + nmps, wi :: wis, wi +. zz))
      ([], [], 0, [], 0.) (rev nonmut_us) in
    let piis = pi0 :: (map (fun wi -> wi /. zz *. (1. -. pi0)) wis) in
    let counts = mut_ct :: counts in
    (tot +. Gsl.Randist.multinomial_lnpdf ~p:(aofl piis) ~n:(aofl counts),
     (gen, mut_ct + nmps, mutty :: types, counts, piis) :: ress) in
  let (ll, ress) = fold sum_gen (0., []) annupdates in
  (ll, rev ress)

let likelihood_diffll ss indty mutp annuds =
  (* let zeros = Array.make (Array.length ss) 0. in *)
  let (obs, exs, ress) = diffll_ress ss ss indty mutp annuds in
  let ll = Gsl.Randist.multinomial_lnpdf ~p:exs ~n:obs in
  (ll, ress)

(* Warning: Censored timeseries should never explicitly list counts such as 0
   for censored counts
  (* lb imputes missing zeros in timeseries *)
  let lb_auds = annupdate_data_nf Mu.cons Mu.cons timeseries in
  (* ub imputes missing cencounts in timeseries *)
  let ub_auds = annupdate_data_nf Mu.cons
    (fun (ty, (fc, tc)) tss -> (ty, (fc, cenct)) :: tss) timeseries in 
  let cenfreq = max_censored_freq cenct (pop_sizes_uds uds) in
  let mutp nnt (ty, (frc, toc)) = norm frc nnt < cenfreq in
  dd is total number of ss including the wildtype s. 
    *)
let ml_mm mutty mutp indty dd annupdates =
  let zeros = Array.make dd 0. in

  let nuest = ml_numigration_rate mutp annupdates in

  let mm sms =
    let (obs, exs, ress) = diffll_ress zeros sms indty mutp annupdates in
    let ll = Gsl.Randist.multinomial_lnpdf 
      ~p:(Array.map2 (fun ss ctz -> if is_nan ss then 0. else exp ss *. ctz)
         sms exs) 
      ~n:obs in
    let optss = Array.map2 (fun ob ex -> log (fl ob /. ex)) obs exs in
    let smean = Mu.meanf_n (atol optss) in
    let optss = Array.map (fun x -> x -. smean) optss in
    (optss, ll, ress) in

  let rec fix fn this_ss last_ss last_ll last_diffress = 
    (* ll and diffress refer to this, next should do a better job *)
    let (next, this_ll, this_diffress) = fn this_ss in
    if is_nan this_ll then failwith "likelihood is nan" ;
    (* let wty = w_pwcty indty next in
    let (lhd, lhdress) = likelihood_cen wty mutty mutp annupdates in *)
    if eq_sels_a next this_ss || this_ll <= last_ll
    then (last_ss, last_ll, last_diffress)
    else fix fn next this_ss this_ll this_diffress in

  let (sels, ml_ll, dre) = fix mm zeros zeros neg_infinity [] in
  (nuest, sels, ml_ll, dre)

let fprint_residuals ouch pre print_ty ts_residuals =
  Printf.fprintf ouch 
    "%sgen\tinit_ps\tfin_ps\ttype\tinitf\tobsf\texpf\tnexpf\trs_inc\n" pre ;
  ignore (fold (fun (initps, finps, (* tys, *) initfs, obsfs, expfs, nexpfs) ngen ->
      Mu.iter4 (fun (* ty *) initf obsf expf nexpf ->
          Printf.fprintf ouch "%s%d\t%d\t%d\tpercents\t%F\t%F\t%F\t%F\t%F\n" pre
            ngen initps finps (* (print_ty ty) *) initf obsf expf nexpf
            ((obsf -. expf) /. sqrt(expf *. (1. -. expf))))
        initfs obsfs expfs nexpfs ;
      ngen + 1)
    1 ts_residuals )

let fprint_ts_residuals ouch pre print_ty ts_residuals =
  Printf.fprintf ouch 
    "%sgen\tinit_ps\tfin_ps\ttype\tinitf\tobsf\texpf\tnexpf\trs_inc\n" pre ;
  ignore (iter 
    (fun (gen, initps, finps, tys, initfs, obsfs, expfs, nexpfs) ->
      Mu.iter5 (fun ty initf obsf expf nexpf ->
          Printf.fprintf ouch "%s%d\t%d\t%d\t%s\t%F\t%F\t%F\t%F\t%F\n"
            pre gen initps finps (print_ty ty) initf obsf expf nexpf 
            ((obsf -. expf) /. sqrt(expf *. (1. -. expf))))
        tys initfs obsfs expfs nexpfs)
    ts_residuals)

let print_residuals tag ts_residuals =
  fprint_residuals stdout (Printf.sprintf "%s\t" tag) ts_residuals

let mutant_count_dist updates =
  let kons ups ns =
    let kons (fc, tc) ns = if fc == 0 then tc :: ns else ns in
    fold kons ns ups in
  fold kons [] updates

let print_mutant_count_dist tag counts =
  Printf.printf "%s\tcount\n%!" tag ;
  iter (fun ct -> Printf.printf "%s\t%d\n%!" tag ct) counts

(* }}} *)

let parse_tsv_raw path =
  let inch = open_in path in
  let rec parse_lines lines linenum ncols = 
    let this_line = try Some (input_line inch) with End_of_file -> None in
    match this_line with
        Some this_line ->
          let vals = Mu.splittab this_line in
          let vlen = List.length vals in
          if vlen <> ncols then failwith 
            (Printf.sprintf "%s:%d: %d cols expected, found %d" 
              path linenum ncols vlen) else
          parse_lines (vals :: lines) (linenum + 1) ncols
      | None -> lines in
  let first_line = try input_line inch 
    with End_of_file -> failwith
      (Printf.sprintf "%s contains no header" path) in
  let header = Mu.splittab first_line in
  rev (parse_lines [header] 2 (List.length header))

let parse_nlsv_ch interp inch = 
  let rec parse_lines vals linenum =
    let this_line = try Some (input_line inch) with End_of_file -> None in
    match this_line with
        Some this_line -> 
          let this_val = try interp this_line 
          with e -> failwith (Printf.sprintf 
            "Failure line %d: %s" linenum (Printexc.to_string e)) in
          parse_lines (this_val :: vals) (linenum + 1)
      | None -> vals in
  rev (parse_lines [] 1)

let parse_nlsv interp path = 
  parse_nlsv_ch interp (open_in path)

let gen_type_count_to_timeseries tuples =
  let kons (gen,typ,cou) mmap = Mp.addf gen (Mp.add typ cou) Mp.empty mmap in
  let mapmap = fold kons Mp.empty tuples in
  let ff (k, tcm) = (k, clean (Mp.bindings tcm)) in
  map ff (Mp.bindings mapmap) (* bindings comes out sorted by (fst, ..)*)

let births_gtc_to_timeseries agef gtcs =
  let kons (gen,typ,cou) mmap = Mp.addf gen (Mp.add typ cou) Mp.empty mmap in
  let gtcsmmb = fold kons Mp.empty gtcs in
  let gens = aofl (map fst (Mp.bindings gtcsmmb)) in
  let ng = Array.length gens in
  let gtcsmm = fold (fun ind gtcsmm ->
      let births = Mp.find gens.(ind) gtcsmmb in
      fold (fun (typ, cou) gtcsmm ->
          let lifetime = agef () in
          snd (Mu.rec_n lifetime (fun (ind, gtcsmm) ->
              if ind < ng then
                (ind + 1,
                 Mp.addf gens.(ind) 
                   (fun mp -> Mp.addf typ ((+) cou) 0 mp) Mp.empty gtcsmm)
              else (ind, gtcsmm))
            (ind, gtcsmm)))
        gtcsmm (Mp.bindings births))
    Mp.empty (Mu.range 0 (ng - 1)) in
  let ff (gen, tcm) = (gen, clean (Mp.bindings tcm)) in
  Mu.drop (agef () - 1) (map ff (Mp.bindings gtcsmm))

let parse_timeseries_gtc path =
  let (header, rows) = match parse_tsv_raw path with
      hd :: rw :: rws -> (hd, rw :: rws)
    | _ -> failwith (Printf.sprintf "%s: no data" path) in
  let header = map (function
      "gen" -> `Gen
    | "type" -> `Type
    | "count" -> `Count
    | _ -> `Other) header in
  let kons2 head value (gen,typ,cou) = match head with
      `Gen -> (int_of_string value, typ, cou)
    | `Type -> (gen, value, cou)
    | `Count -> (gen, typ, int_of_string value)
    | _ -> (gen,typ,cou) in
  map (fun row -> Mu.fold2 kons2 (0, "", 0) header row) rows

let parse_timeseries path =
  gen_type_count_to_timeseries (parse_timeseries_gtc path)

let parse_pss path = parse_nlsv int_of_string path

let parse_tsv path years =
  let yeardata year =
    let ic = Scanf.Scanning.open_in
      (Printf.sprintf "%s%d.tsv" path year) in
    let rec scan acc =
      try scan (Scanf.bscanf ic "%s@\t%d\n" (fun a b c -> (a, b) :: c) acc)
      with End_of_file -> acc
        | e -> Printf.eprintf "Error in %d:\n%!" year ; raise e in
    rev (scan []) in
  map yeardata years

let parse_tsv_date_range path =
  let dates = Mu.slurp_stdout
    (Printf.sprintf "ls %s*.tsv | sed -e 's/.tsv//g' | sed -e 's!%s!!g'"
      path path) in
  let dates = map int_of_string (Mu.splitnl dates) in
  Mu.min_max dates

type params = {
  mu : (float * float) ;
  ne : (float * float) ;
  sels : (float array * float array) }


let fprint_bootstrap ouch bootstrap =
  Printf.fprintf ouch "param\trun\tind\tval\n" ;
  ignore (fold (fun (mu, sels, srep, ne, _) run ->
      Printf.fprintf ouch "mu\t%d\tNA\t%F\n" run mu ;
      Printf.fprintf ouch "srep\t%d\tNA\t%F\n" run srep ;
      Printf.fprintf ouch "ne\t%d\tNA\t%F\n" run ne ;
      Array.iteri (fun ind sel ->
          Printf.fprintf ouch "s%d\t%d\t%d\t%F\n" (ind + 1) run ind sel) sels ;
      run + 1) 0 bootstrap)

let parse_params path =
  let (header, rows) = match parse_tsv_raw path with
      hd :: rw :: rws -> (hd, rw :: rws)
    | _ -> failwith (Printf.sprintf "%s: no data" path) in
  let konscell header x (par, ind, vvv, var) = match header with
      "param" -> (Some x, ind, vvv, var)
    | "ind" ->   (par, Some x, vvv, var)
    | "val" ->   (par, ind, Some x, var)
    | "var" ->   (par, ind, vvv, Some x)
    | _ ->       (par, ind, vvv, var) in
  let params = map (fun row -> Mu.fold2 konscell (None, None, None, None) 
                                 header row)
                 rows in
  let sparams = Mu.filtmap (fun (_, ind, vvv, var) ->
      try Some (int_of_string (Mu.opt_req ind), vvv, var) 
      with 
        Mu.Required_arg -> failwith "ind column missing"
      | Failure _ -> None) params in
  let dd = List.length sparams in 
  let (mu, mu_var) = match Mu.filtmap (fun (par, _, vvv, var) -> 
        if par = Some "mu" then Some (vvv, var) else None) params with
      [(Some mu, Some mu_var)] -> (float_of_string mu, float_of_string mu_var)
    | _ -> failwith (Printf.sprintf "No parameter mu in %s" path) in
  let ne = match Mu.filtmap (fun (par, _, vvv, var) -> match par, vvv, var with
            Some "ne", Some vvv, Some var -> 
              Some (vvv, var)
          | _ -> None) params with
      [(a, b)] -> let fs = float_of_string in (fs a, fs b)
    | _ -> failwith (Printf.sprintf "No parameter ne in %s" path) in
  let (ss, svars) = (Array.make dd 0., Array.make dd 0.) in
  iter (fun (ind, vvv, var) -> match (vvv, var) with 
          (Some vv, Some vr) -> 
            ss.(ind) <- float_of_string vv ; 
            svars.(ind) <- float_of_string vr
        | _ -> failwith (Printf.sprintf "val or var absent for index %d" ind))
    sparams ;
  { mu=(mu, mu_var); ne=ne; sels = (ss, svars) }

let params_of_mu_sels indty mu sels annuds lhdress = 
  let (var_mu, var_ss) = observed_information_var mu indty sels annuds in
  let ne = ne_estimate (obsexps_of_lhdress lhdress) in
  { mu = (mu, var_mu) ; ne = ne ; sels = (sels, var_ss) }

(* Resample transitions each generation so as to preserve the source
  population size every generation under expectation. This means that,
  typically,
    pop_sizes_uds ups != pop_sizes_uds (resample ups)
  although the two are equal under expectation.  Furthermore the resampled
  updates are somewhat incoherent from generation to generation; in particular
  pop_sizes_uds rups and pop_sizzes_uds_prime rups will typically bear no
  resemblence to one another although the two are equal under expectation as
  well. *)
let resample annuds =
  let resamp (gen, ups) =
    let upsa = aofl ups in
    let nups = Array.length upsa and ps = pop_size_a ups in
    let maxfst = Mu.max (map (fun (ty, (fc, tc)) -> fc) ups) in
    (* Are we done drawing transitions? *)
    let done_drawing_p cur =
      let dd = ps - cur in (* gap left to fill *)
      if maxfst <= dd then false (* always keep drawing in this case *)
      else let exp = norm (Mu.sum (* expectation after one more draw *)
        (map (fun (_,(ct,_)) -> if ct <= dd then ps else cur + ct) ups)) 
             nups in
      let probf = fl (ps - cur) /. (exp -. fl cur) in
      Gsl.Randist.bernoulli Rand.rng ~p:probf == 0 in
    let unifa = Array.make nups (1. /. fl nups) in
    let unif = Gsl.Randist.discrete_preproc unifa in (* uniform distribution *)
    (gen, List.sort compare 
      (snd (Mu.rec_p
        (fun (totct, _) -> done_drawing_p totct) (* when to quit *)
        (fun (totct, acc) -> (* draw another *)
          let (_,(fc,tc)) as up = upsa.(Gsl.Randist.discrete Rand.rng unif) in
          (totct + fc, up :: acc))
        (0, [])))) in
  map resamp annuds

(* let bootstrap_params_pwc nruns infer aupds =
  Printf.printf "Bootstrap resample %!" ;
  let res = Mu.rec_n nruns (fun knil ->
    Printf.printf ".%!" ;
    let aupds = resample aupds in
    let (mu, sels, srep, ne, ll, _, _) = infer aupds in
    (mu, sels, srep, ne, ll) :: knil) [] in
  Printf.printf "\n%!" ; res *)

let bootstrap_params_pwc_inc_fprint ch nruns infer aupds =
  Printf.printf "Bootstrap resample %!" ;
  Printf.fprintf ch "param\trun\tind\tval\n" ;
  let res = Mu.rec_n nruns (fun (run, knil) ->
    Printf.printf ".%!" ;
    let aupds = resample aupds in
    let (mu, sels, srep, ne, ll, _, _) = infer aupds in
    Printf.fprintf ch "mu\t%d\tNA\t%F\n" run mu ;
    Printf.fprintf ch "srep\t%d\tNA\t%F\n" run srep ;
    Printf.fprintf ch "ne\t%d\tNA\t%F\n" run ne ;
    Array.iteri (fun ind sel ->
        Printf.fprintf ch "s%d\t%d\t%d\t%F\n" (ind + 1) run ind sel) sels ;
    Printf.fprintf ch "%!" ;
    (run + 1, (mu, sels, srep, ne, ll) :: knil)) (0, []) in
  Printf.printf "\n%!" ; res

let infer_params wtty mutty mutp indty nbins annuds =
  let (mu, sels, ll, dre) = 
    ml_mm mutty mutp indty nbins annuds in
  let wtann = match wtty with
    Some wt -> aggregate_wt wt indty annuds | None -> annuds in
  let (_, lre) = likelihood_cen mu (w_pwcty indty sels) mutty mutp wtann in
  (params_of_mu_sels indty mu sels annuds lre, ll, dre, lre)

(* infer but don't compute confidence intervals *)
let infer_nocis mutty mutp indty nbins annuds = 
    let (mu, sels, ll, dre) =
      ml_mm mutty mutp indty nbins annuds in
    let (_, lre) = likelihood_cen mu (w_pwcty indty sels) mutty mutp annuds in
    let ne = fst (ne_estimate (obsexps_of_lhdress lre)) in
    let srep = num_replacement_fitness mu indty sels mutp annuds in
    (mu, sels, srep, ne, ll, dre, lre)

let fprint_params ouch params minfreq wtty mutty mutp indty annuds =
  let array_tl arr = Array.sub arr 1 (Array.length arr - 1) in
  let (mu, mu_var) = params.mu in
  let (ne, ne_var) = params.ne in
  let (sels, sels_var) = params.sels in
  let (lcis, ucis) = vars_to_cis sels sels_var in
  let (sis, sis_var) = (match wtty with None -> (sels, sels_var) 
    | Some _ -> (array_tl sels, array_tl sels_var)) in
  let (silcis, siucis) = (match wtty with None -> (lcis, ucis) 
    | Some _ -> (array_tl lcis, array_tl ucis)) in
  let srep = num_replacement_fitness mu indty sels mutp annuds in
  let wtannuds = match wtty with
    Some wt -> aggregate_wt wt indty annuds | None -> annuds in
  let wty_fdsel = w_pwcty indty sels in
  let (ll, lress) = likelihood_cen mu wty_fdsel mutty mutp wtannuds in
  let (s0ll, s0ress) = likelihood_cen mu wty_neutral mutty mutp wtannuds in
  let (s0ne, s0ne_var) = ne_estimate (obsexps_of_lhdress s0ress) in
  let (minf, maxf) = min_max_freq_nowt wtty annuds in
  let mw = 1.96 *. sqrt (mu_var) in
  let neciw = 1.96 *. sqrt (ne_var) in
  let s0neciw = 1.96 *. sqrt (ne_var) in
  Printf.fprintf ouch "param\tind\tval\tvar\tlci\tuci\n%!" ;
  Printf.fprintf ouch "mu\tNA\t%F\t%F\t%F\t%F\n%!" 
    mu mu_var (mu-.mw) (mu+.mw) ;
  Printf.fprintf ouch "ne\tNA\t%F\t%F\t%F\t%F\n%!" 
    ne ne_var (ne-.neciw) (ne+.neciw) ;
  Printf.fprintf ouch "maxf\tNA\t%e\tNA\t%e\t%e\n%!" maxf maxf maxf ;
  ignore (Mu.opt_fmap (fun minf -> 
    Printf.fprintf ouch "minf\tNA\t%e\tNA\t%e\t%e\n%!" minf minf minf) 
    minfreq) ;
  Printf.fprintf ouch "dminf\tNA\t%e\tNA\t%e\t%e\n%!" minf minf minf ;
  Printf.fprintf ouch "srep\tNA\t%F\tNA\t%F\t%F\n%!" srep srep srep ;
  (match wtty with None -> () | Some _ ->
    Printf.fprintf ouch "swt\tNA\t%F\t%F\t%F\t%F\n%!" 
      sels.(0) sels_var.(0) lcis.(0) ucis.(0));
  Array.iteri (fun ind sel ->
      Printf.fprintf ouch "s%d\t%d\t%F\t%F\t%F\t%F\n%!" 
        (ind + 1) ind sel sis_var.(ind) silcis.(ind) siucis.(ind)) sis ;
  Printf.fprintf ouch "ll\tNA\t%F\tNA\tNA\tNA\n%!" ll ;
  Printf.fprintf ouch "s0ll\tNA\t%F\tNA\tNA\tNA\n%!" s0ll ;
  Printf.fprintf ouch "s0ne\tNA\t%F\t%F\t%F\t%F\n%!" 
    s0ne s0ne_var (s0ne-.s0neciw) (s0ne+.s0neciw) ;
  Printf.fprintf ouch "s0lrpval\tNA\t%e\tNA\tNA\tNA\n%!" 
    (Rand.llr_pval (Array.length sels - 1) (ll -. s0ll));
  (match wtty with None -> () | Some wtt ->
    let (indf, breaks) = log_binning 1 wtty minfreq annuds in
    let (totbins, indty) = parse_indty_args 1 indf wtty minfreq in
    let (params,_,_,_) = infer_params wtty mutty mutp indty totbins annuds in
    let (mu,_) = params.mu in
    let (swt0sels,_) = params.sels in
    let wty_swt0 = w_pwcty indty swt0sels in
    let (swt0ll,_) = likelihood_cen mu wty_swt0 mutty mutp wtannuds in (
    Printf.fprintf ouch "swt0ll\tNA\t%F\tNA\tNA\tNA\n%!" swt0ll ;
    Printf.fprintf ouch "swt0lrpval\tNA\t%e\tNA\tNA\tNA\n%!" 
       (Rand.llr_pval (Array.length sels - Array.length swt0sels) 
         (ll -. swt0ll))))
