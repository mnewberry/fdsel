
(* Use GSL's random number generator *)
let rng = Gsl.Rng.make (Gsl.Rng.default ())

let time_of_day () = int_of_float (Unix.gettimeofday () *. 1000000.)

let init seed = Gsl.Rng.set rng (Nativeint.of_int seed)
let int i = if i = 1 then 0 else Gsl.Rng.uniform_int rng i
let float scale = scale *. Gsl.Rng.uniform rng
let lcalphastr n = String.map (fun _ -> Char.chr (int 26 + 97)) 
  (String.make n ' ')

let binom n p = Gsl.Randist.binomial rng p n
let multinom n ps = Array.to_list
  (Gsl.Randist.multinomial rng n (Array.of_list ps))

(* Use ocaml's random number generator
let init = Random.init
let int i = if i = 1 then 0 else Random.int i
let float = Random.float  *)

let exp lambda = log (float 1.0) /. (0. -. lambda) 

(* This gives the cutoff value for the LRT with df degrees of freedom.
   Compare this value to the difference in log-likeihoods. *)
let lrt_cutoff df = Gsl.Cdf.chisq_Qinv 0.05 (float_of_int df) /. 2.
let llr_pval df x = Gsl.Cdf.chisq_Q (2. *. x) (float_of_int df)
