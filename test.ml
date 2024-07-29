module F = Fdsel_lib
let rev = Mu.rev and hd = List.hd

open OUnit2 ;;

let () = ignore (run_test_tt_main ("all" >::: [
  ("test updates" >::: (
    let ts_gtc = F.parse_timeseries_gtc "test/timeseries.tsv" in
    let ts = F.gen_type_count_to_timeseries ts_gtc in
    let map = Mu.map in
    let uds = F.update_data (map snd ts) in
    let annuds = F.annupdate_data ts in
    let annuds_ub = F.annupdates_ub 20 "CENSORED" ts in
    let (indf, breaks) = F.log_binning 10 None None annuds in
    let indty = F.findty indf in
    let sels = Array.make 10 0. in
    let (obs, exs) = F.diffll sels indf sels uds in
    let (obsa, exsa, ress) = F.diffll_ress 
      sels sels indty F.default_mutp annuds in
    let sels1 = Array.mapi (fun ii _ -> log (obs.(ii) /. exs.(ii))) obs in
    let mu = F.ml_mutation_rate uds in
    let lhd = F.likelihood mu (F.w_pwc indf sels1) uds in
    let (lhdc, lhdress) = F.likelihood_cen mu (F.w_pwcty indty sels1) 
      "MUT" F.default_mutp annuds in
    (* removed legacy test
      let (pmu, psels, pll) = F.ml_mm_pwclf None indf sels uds in 
      let (cmu, csels, cll,_) = F.ml_mm "" F.default_mutp indty 10 annuds in *)
    [("annupdates preserves updates" >:: (fun _ -> assert_equal 
      uds (F.bare_updates annuds)));
    ("annupdates_ub preserves popsizes" >:: (fun _ -> assert_equal 
      (F.pop_sizes_auds annuds) (F.pop_sizes_auds annuds_ub)));
    ("annupdates_ub preserves prime popsizes" >:: (fun _ -> assert_equal 
      (F.pop_sizes_auds_prime annuds) (F.pop_sizes_auds_prime annuds_ub)));
    ("pop_sizes methods agree 1" >:: (fun _ -> assert_equal
      (F.pop_sizes_ts ts)
      ((hd (F.pop_sizes_uds uds)) :: F.pop_sizes_uds_prime uds)));
    ("pop_sizes methods agree 2" >:: (fun _ -> assert_equal
      (F.pop_sizes_ts ts)
      (rev (hd (rev (F.pop_sizes_uds_prime uds)) 
            :: rev (F.pop_sizes_uds uds)))));
    ("diffll_ress agrees with diffll obs when mutp = (frc = 0)" >:: 
      (fun _ -> assert_equal obs (Array.map float_of_int obsa)));
    ("diffll_ress agrees with diffll exp when mutp = (frc = 0)" >:: 
      (fun _ -> assert_equal exs exsa));
    ("likelihoods agree" >:: 
      (fun _ -> assert_equal ~cmp:cmp_float ~printer:string_of_float
        lhd lhdc));
    ("mutation and numigration agree when mutp = f == 0" >:: 
      (fun _ -> assert_equal (F.ml_mutation_rate uds) 
        (F.ml_numigration_rate F.default_mutp annuds)))
    (* removed legacy test
    ("ml_mm functions agree under the circumstances" >:: 
     (fun _ -> assert_equal ~printer:(Std.dump) (pmu, psels) (cmu, csels)));
    ("ml_mm functions agree in ll up to numerical error" >:: 
     (fun _ -> assert_equal ~cmp:(cmp_float ~epsilon:1e-15) 
        ~printer:string_of_float pll cll)) *)
        ]));
  ("test explicit_indf" >::: (
    let breaks = [0.1;0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9] in
    ("lower-bound" >:: 
      (fun _ -> assert_equal 0 (F.explicit_indf breaks 0.0))) :: 
    ("boundary-value" >::
      (fun _ -> assert_equal 0 (F.explicit_indf breaks 0.1))) ::
    ("upper-bound" >::
      (fun _ -> assert_equal 9 (F.explicit_indf breaks 1.0))) ::
    (snd (Mu.rec_n 8 (fun (i, kn) -> (i + 1, 
       ((Printf.sprintf "interim test %d" i) >::
         (fun _ -> assert_equal i (F.explicit_indf breaks 
           (float_of_int i /. 10. +. 0.01)))) :: kn)) (1, [])))));
  ("indf and bin counts" >::: (
    let nbins = 10 in
    let breaks = [0.1;0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9] in
    let sels = [|0.1;0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9;1.0|] in
    let test_indf indf x sel = assert_equal ~printer:string_of_float 
      (exp sel) (F.w_pwc indf sels x) in [
    ("equal_bins_ind" >::: [
      ("0." >::
        (fun _ -> test_indf (F.equal_bins_ind None nbins) 0.0001 0.1));
      ("1." >::
        (fun _ -> test_indf (F.equal_bins_ind None nbins) 1. 1.0));
      ("too high" >:: (fun _ -> assert_raises 
        (Failure "Frequency outside range")
        (fun () -> F.equal_bins_ind None nbins 5.)));
      ("too low" >:: (fun _ -> assert_raises
        (Failure "Frequency outside range")
        (fun () -> F.equal_bins_ind None nbins (-1.))))]);
    ("log_bins" >::: [
      ("0." >::
        (fun _ -> test_indf (fst (F.log_bins nbins 0.1 0.9)) 0.0001 0.1));
      ("1." >::
        (fun _ -> test_indf (fst (F.log_bins nbins 0.1 0.9)) 1. 1.0));
      ("too high" >::
        (fun _ -> test_indf (fst (F.log_bins nbins 0.1 0.9)) 5. 1.0));
      ("too low" >::
        (fun _ -> test_indf (fst (F.log_bins nbins 0.1 0.9)) (-1.) 0.1))]);
    ("explicit_bins" >::: [
      ("0." >::
        (fun _ -> test_indf (F.explicit_indf breaks) 0.0001 0.1));
      ("1." >::
        (fun _ -> test_indf (F.explicit_indf breaks) 1. 1.0));
      ("too high" >::
        (fun _ -> test_indf (F.explicit_indf breaks) 5. 1.0));
      ("too low" >::
        (fun _ -> test_indf (F.explicit_indf breaks) (-1.) 0.1))])]));
  ("gen_intervals" >::: (
    let ts1 = F.gen_type_count_to_timeseries 
      (F.parse_timeseries_gtc "test/g1ts.tsv") in
    let ts2 = F.gen_type_count_to_timeseries 
      (F.parse_timeseries_gtc "test/g12ts.tsv") in
    [("single generation" >:: 
      (fun _ -> assert_equal [] (F.gen_intervals [List.hd ts1])));
    ("no generation" >:: (fun _ -> assert_equal [] (F.gen_intervals [])));
    ("generation 1" >:: (fun _ -> assert_equal [1] (F.gen_intervals ts1)));
    ("odd generations" >::
      (fun _ -> assert_equal [1;2] (F.gen_intervals ts2)))]))
        ]))
