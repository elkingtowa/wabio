;;; Jeff Shrager's General Lisp Numerical Utilities
;;; Copyright (c) Jeff Shrager 1999 2000 2001 
;;; Contact: jshrager@andrew2.stanford.edu; 650.325.1521.287

; (load "g:/jshrager/lib/stats.lsp")
; (compile-file "g:/jshrager/lib/stats.lsp")

(defmacro display (&rest l)
  `(progn ,@(loop for i in l
		  collect `(format t "~a = ~a~%" ',i ,i)))) 


;;; --- Various math and stats routines.

;;; Returns a string that is rounded to the appropriate number of
;;; digits, but the only thing you can do with it is print it.  It's
;;; just a convenience hack for rounding recursive lists.

(defun pround (n v)
  (if (listp v)
      (mapcar #'(lambda (v) (pround n v)) v)
    (if (numberp v)
	(format nil (format nil "~~,~af" n) v)
      v)))

;;; And some conveniences:

(defun p2 (v) (pround 2 v))
(defun p1 (v) (pround 2 v))
(defun p3 (v) (pround 3 v))
(defun p2* (l) (mapcar #'p2 l))

;;; One way t-test to see if a group differs from a numerical mean
;;; target value.  From Misanin & Hinderliter p. 248.

(defun t1-test (values target &optional (warn? t))
  (let ((t1-value (t1-value values target)))
    (values t1-value
	    (t-p-value (abs t1-value) (1- (length values)) warn?))))

(defun t1-value (values target)
  (/ (- (mean values) target)
     (/ (standard-deviation values)
	(sqrt (length values)))))

;;; Select n random sublists from a list, without replacement.  This
;;; copies the list and then destroys the copy.  N better be less than
;;; or equal to (length l).

(defun n-random (n l &aux r)
  (setq l (copy-list l))
  (dotimes (k n)
    (if (null (cdr l)) ; this has to be the last one
        (progn (push (car l) r)
	       (return nil))
      ;; If we're taking the first, it's easy, otherwise 
      ;; the first item get shoved into the middle.
      (let ((take (random (length l))))
        (if (zerop take)
	    (push (car l) r)
	  (let ((nc (nthcdr take l)))
	    (push (car nc) r)
	    (rplaca nc (car l))
	    )
	  )
	(pop l)
	) ; let
      ) ; if 
    ) ; dotimes
  r)

;;; T2-test calculates an UNPAIRED t-test.
;;; From Misanin & Hinderliter p. 268.   The t-cdf part is inherent in 
;;; xlispstat, and I'm not entirely sure that it's really the right
;;; computation since it doens't agree entirely with Table 5 of M&H, but
;;; it's close, so I assume that M&H have round-off error.

(defun t2-test (l1 l2)
  (if (or (equal l1 l2) ; protect against 0 result from sdiff!
	  (null (cdr l1)) ; and other stupidity!
	  (null (cdr l2)))
    (values "n/a" "n/a")
    (let ((tval (t2-value l1 l2)))
      (values tval (t-p-value (abs tval) (- (+ (length l1) (length l2)) 2))))))

(defun t2-value (l1 l2)
  (let* ((n1 (float (length l1)))
	 (n2 (float (length l2)))
	 (s21 (s2 l1 n1))
	 (s22 (s2 l2 n2))
	 (m1 (mean l1))
	 (m2 (mean l2))
 	 (s2p (/ (+ (* s21 (1- n1))
		    (* s22 (1- n2)))
		 (+ (1- n1) (1- n2))))
	 )
    (/ (- m1 m2)
       (sqrt (+ (/ s2p n1) (/ s2p n2))))
    ))

(defun s2 (l n)
  (/ (- (sum (mapcar #'(lambda (a) (* a a)) l))
	(/ (let ((r (sum l))) (* r r)) n))
     (1- n)))

;;; For non-xlispstat, we compute t-cdf the hard way.  Sometime I
;;; ought to put in a real computation, but for the moment, it's a
;;; total cheat, just returning the symbol '>.05 or '<.05.  Each
;;; Entry is a df followed by the critical point at that df.

(defvar *t-cdf-critical-points-table-for-.05*
  '((1 . 6.31)
    (2 . 2.92)
    (3 . 2.35)
    (4 . 2.13)
    (5 . 2.02)
    (6 . 1.94)
    (7 . 1.89)
    (8 . 1.86)
    (9 . 1.83)
   (10 . 1.81)
   (15 . 1.75)
   (20 . 1.72)
   (25 . 1.71)
   (30 . 1.70)
   (40 . 1.68)
   (60 . 1.67)
  (120 . 1.66)))

(defun t-p-value (x df &optional (warn? t))
  (if (> df 120)
      (progn 
	(if warn? 
	    (format t "Warning df (~a) is off the scale.  Using df=120~%" df))
	(setq df 120)))
  (dolist (tcp *t-cdf-critical-points-table-for-.05*)
    (if (or (= df (car tcp))
	    (< df (car tcp)))
	;;                                      nil???
	(return (if (> x (cdr tcp)) '*/p<.05! 'ns/p>.05))))
  )

;;; --- This emulates some of the xlispstat functions, and other stats
;;; utilities. 

;;; Some math ops that take many args.

(defun sum (l &aux (sum 0))
  (dolist (n l) (incf sum n)) sum)

(defun sqr (a)
  (if (numberp a) (expt a 2)
      (mapcar #'* a a)))

;;; This version of min and max operate on a list if there's only one arg.

(defun max* (l &rest ll &aux m)
  (if ll (setq l (cons l ll)))
  (setq m (pop l))
  (dolist (i l)
    (if (> i m) (setq m i)))
  m)

(defun min* (l &rest ll &aux m)
  (if ll (setq l (cons l ll)))
  (setq m (pop l))
  (dolist (i l)
    (if (< i m) (setq m i)))
  m)

;;; Warning, you can't use apply on most of the math operations
;;; because the arglist limits a too small for any large list.

(defun mean (l)
  (/ (sum l) (float (length l))))

;;; Computes a mean protected where there will be a divide by zero, and gives us n/a in that case.

(defun protected-mean (l)
  (loop with sum = 0
	with n = 0
	as value in l
	when (numberp value)
	do (incf sum value) (incf n)
	finally (return (cond ((not (zerop n)) (/ sum (float n)))
			      (t 'n/a)))))

;;; This is specially protected for zero divisors.

(defun standard-deviation (l)
  (if (null (cdr l))
      0.0
    (sqrt
     (let ((m (mean l)))
       (* (/ 1.0 (1- (length l)))
	  (sum (mapcar #'(lambda (x) (expt (- x m) 2)) l))
	  )))))

(defun standard-error (l)
  (/ (standard-deviation l) (sqrt (length l))))

;;; Lmean takes the mean of entries in a list of lists vertically.  So:
;;;  (lmean '((1 2) (5 6))) -> (3 4)  The args have to be the same length.

(defun lmean (ll)
  (let* ((n (length ll))
	 (sums (copy-list (pop ll))) ; copy so that incf doesn't wreck us!
	 (len (length sums)))
    (dolist (l ll)
      (dotimes (k len)
	(incf (nth k sums) (nth k l))
	))
    (setq n (float n))
    (mapcar #'(lambda (v) (/ v n)) sums)
    )
  )

;;; --- Two-Way Anova.  (From Misanin & Hinderliter, 1991, p. 367-) This
;;; is specialized for four groups of equal n, called by their plot
;;; location names: left1 left2 right1 right2.

(defun anova2 (a1b1 a1b2 a2b1 a2b2)
  (let* ((n (length a1b1)) ; All are the same, so any will do here.
         (npq (* 4 n)) ; This is specialized for a 2x2 design; always 4 levels.
	 (a1 (append a1b1 a1b2))
	 (suma1 (sum a1))
	 (a2 (append a2b1 a2b2))
	 (suma2 (sum a2))
	 (b1 (append a1b1 a2b1))
	 (sumb1 (sum b1))
	 (b2 (append a1b2 a2b2))
	 (sumb2 (sum b2))
	 (allscores (append a1 a2))
	 (sym1 (float (/ (sqr (sum allscores)) npq)))
	 (sym2 (float (sum (sqr allscores))))
	 (np (* n 2))
	 (sym3 (float (/ (+ (sqr suma1) (sqr suma2)) np)))
	 (nq np)
	 (sym4 (float (/ (+ (sqr sumb1) (sqr sumb2)) np)))
	 (sym5 (float (/ (+ (sqr (sum a1b1)) (sqr (sum a2b1)) 
			    (sqr (sum a1b2)) (sqr (sum a2b2)))
			 n)
		      ))
	 (ssbg (- sym5 sym1))
	 (ssa (- sym3 sym1))
	 (ssb (- sym4 sym1))
	 (ssab (+ sym1 (- (- sym5 sym4) sym3)))
	 (sswg (- sym2 sym5))
	 (sstot (- sym2 sym1))
	 (df1 3)
	 (df2 1)
	 (df3 1)
	 (df4 1)
	 (df5 (* 4 (1- n)))
	 (dftot (1- npq))
	 (msbg (float (/ ssbg df1)))
	 (msa (float (/ ssa df2)))
	 (msb (float (/ ssb df3)))
	 (msab (float (/ ssab df4)))
	 (mswg (float (/ sswg df5)))
	 )
    (list :ssbg ssbg :dfbg df1 :msbg msbg :fbg (/ msbg mswg)
	  :ssa ssa :dfa df2 :msa msa :fa (/ msa mswg)
	  :ssb ssb :dfb df3 :msb msb :fb (/ msb mswg)
	  :ssab ssab :dfab df4 :msab msab :fab (/ msab mswg)
	  :sswg sswg :dfwg df5 :mswg mswg
	  :sstot sstot :dftot dftot
	  :a1b1 a1b1 :a1b2 a1b2 :a2b1 a2b1 :a2b2 a2b2)
    ))

(defun all-squares (as bs &aux squares)
  (dolist (a as)
    (dolist (b bs)
      (push (sqr (* a b)) squares)
      )))

(defun testanova2 ()
  (let ((a1b1 '(1 2 3 4 5))
	(a1b2 '(3 4 5 6 7))
	(a2b1 '(5 6 7 8 9))
	(a2b2 '(4 5 6 7 8))
	)
    (anova2 a1b1 a1b2 a2b1 a2b2)))

#| From VassarStats:

(In the above, A are rows, B are cols)

--------------------------------------------------------------------------------
ANOVA SUMMARY
Source   SS     df       MS     F        P 
bg       43.75   3   
rows     31.25   1    31.25  12.5 0.002749 
columns   1.25   1     1.25   0.5 0.489672 
r x c    11.25   1    11.25   4.5 0.049865 
wg       40     16     2.5  
Total    83.75  19   

(bg = between groups; wg = within groups (error))


Results are summarised in ANOVA summary table:

Source Sum of Squares Df Mean Squares  F-ratio p 
Age 11.25 1 11.25 37.5 .000 
Sex 101.25 1 101.25 337.5 .000 
Age * sex 1.25 1 1.25 4.167 .058 
Error 4.8 16 .3     
Total 118.55 20       

1. SS for factors and total calculated by formula
2. SS (error) = SS (total) – SS (factors): 118.55 – 11.25 – 101.25 – 1.25 = 4.8
Df calculated from N and number of conditions 
MS calculated: SS / DF (for each factor) 
Age: 11.25 / 1 = 11.25
Sex: 101.25 / 1 = 101.25
Interaction Age * Sex (not notation): 1.25 /1 = 1.25
Error: 4.8 / 16 = 0.3
F ratio calculated: MS (factor) / MS (error) 
Age: 11.25 / 0.3 = 37.5
Sex: 101.25 / 0.3 = 337.5
Interaction Age * Sex: 1.25 / 0.3 = 4.167

How to report:

Main effect of Age: F(1,16) = 37.5, p < 0.0001
Main effect of Sex: F(1,16) = 337.5, p < 0.0001
Interaction between Age and Sex (or: Age x Sex interaction): F(1,16) = 4.167, p = 0.058



|#

;;; --- Two way ANOVA with repeated measures on one dimension.  From
;;; Ferguson & Takane, 1989, p. 359.  Data is organized differently
;;; for this test.  Each group (g1 g2) contains list of all subjects'
;;; repeated measures, and same for B.  So, A: ((t1s1g1 t2s1g1 ...)
;;; (t1s2g2 t2s2g2 ...) ...)  Have to have the same number of test
;;; repeats for each subject, and this assumes the same number of
;;; subject in each group.

(defun anova2r (g1 g2)
  (let* ((c (length (car g1)))
	 (r 2) ; only two groups in this special case
	 (tsr1 (mapcar #'sum g1))
	 (t1 (sum tsr1))
	 (t1c (let (temp)
		(dotimes (n c temp)
		 (push (sum (mapcar #'(lambda (s) (nth n s)) g1)) temp))))
	 (tsr2 (mapcar #'sum g2))
	 (t2 (sum tsr2))
	 (t2c (let (temp)
		(dotimes (n c temp)
		 (push (sum (mapcar #'(lambda (s) (nth n s)) g2)) temp))))
	 (tc (mapcar #'+ t1c t2c))
	 (total (+ t1 t2))
	 (n (length g1))
	 (q1 (* (/ 1.0 c) (sum (append (sqr tsr1) (sqr tsr2)))))
	 (q2 (* (/ 1.0 (* c n)) (+ (sqr t1) (sqr t2))))
	 (q3 (* (/ 1.0 (* 2 n)) (sum (sqr tc))))
	 (q4 (* (/ 1.0 n) (sum (append (sqr t1c) (sqr t2c)))))
	 (q5 (sum (append (mapcar #'sum (mapcar #'sqr g1))
			  (mapcar #'sum (mapcar #'sqr g2)))))
	 (q6 (/ (sqr total) (* n 2.0 c)))
	 (ssbs (- q1 q6))
	 (ssr (- q2 q6))
	 (sssr (- q1 q2))
	 (ssws (- q5 q1))
	 (ssc (- q3 q6))
	 (ssrc (+ q6 (- (- q4 q2) q3)))
	 (sscsr (+ q2 (- (- q5 q1) q4)))
	 (sstotal (- q5 q6))
	 (dfr (1- r))
	 (dfc (1- c))
	 (dfrc (* (1- r) (1- c)))
	 (dfsr (* r (1- n)))
	 (dfcsr (* r (1- n) (1- c)))
	 (msr (/ ssr dfr))
	 (mssr (/ sssr dfsr))
	 (msc (/ ssc dfc))
	 (msrc (/ ssrc dfrc))
	 (mscsr (/ sscsr dfcsr))
	 (dftotal (+ dfr dfc dfrc dfsr dfcsr))
	 (fr (/ msr mssr))
	 (fc (/ msc mscsr))
	 (frc (/ ssrc mscsr))
	 )
    (if (not (= (length g1) (length g2)))
	(format t "Warning, ANOVA2R design has unbalanced cells!~%"))
    (list :ssbs ssbs :ssr ssr :sssr sssr :ssws ssws :ssc ssc 
	  :ssrc ssrc :sscsr sscsr :sstotal sstotal
	  :dfr dfr :dfsr dfsr :dfc dfc :dfrc dfrc :dfcsr dfcsr :dftotal dftotal
	  :msr msr :mssr mssr :msc msc :msrc msrc :mscsr mscsr
	  :fr fr :fc fc :frc frc)
    ) ; close let*
  )

(setq an2rd1 '( (2 7 6 7 9)
      (4 3 7 12 14)
      (7 6 4 12 10)
      (1 3 3 6 6)))
(setq an2rd2 '( (4 4 7 9 1)
      (10 12 12 12 16)
      (8 7 8 12 10)
      (5 7 6 7 8)))

;;; One way simple ANOVA, from Neter, et al. p677+.
;;; Data is give as a list of lists, each one representing a treatment, and each containing
;;; the observations.

;;; Example from Neter p.677

(setq neter677data 
      '((11 17 16 14 15) (12 10 15 19 11) (23 20 18 17) (27 33 22 26 28)))

(defun anova1 (d) ; Note that dots are replaced by dashes, so: Y.. = Y-- here
  (let* ((meanYi- (mapcar #'mean d))
	 (serrYi- (mapcar #'standard-error d))
	 (r (length d))
	 (Yi- (mapcar #'sum d))
	 (ni (mapcar #'length d))
	 (Y-- (sum Yi-))
	 (meanY-- (mean meanYi-))
	 (ntotal (sum ni))
	 (SSTO (sum (mapcar #'(lambda (treatment) 
				(sum (mapcar #'(lambda (Oij) 
						 (expt (- Oij meanY--) 2)) treatment))) d)))
	 (SSE (sum (mapcar #'(lambda (treatment meanYi-) 
			       (sum (mapcar #'(lambda (Oij) 
						(expt (- Oij meanYi-) 2)) treatment))) d meanYi-)))
	 (SSTR (sum (mapcar #'(lambda (meanYi- ni) (* ni (expt (- meanYi- meanY--) 2))) meanYi- ni)))
	 (SSTRdf (- r 1))
	 (SSEdf (- ntotal r))
	 (SSTOdf (- ntotal 1))
	 (MSTR (/ SSTR SSTRdf))
	 (MSE (/ SSE SSedf))
	 (F* (/ MSTR MSE))
	 )
	 (list :SUMMARY (format nil "F(.95,~a,~a) = ~a" SSTRdf SSEdf F*) 
	       :r r
	       :meanYi- meanYi-
	       :serrYi- serrYi-
	       :Yi- Yi- 
	       :ni ni
	       :Y-- Y--
	       :meanY-- meanY--
	       :ntotal ntotal
	       :SSTO SSTO
	       :SSE SSE
	       :SSTR SSTR 
	       :SSTRdf SSTRdf
	       :SSEdf SSEdf
	       :SSTOdf SSTOdf
	       :MSTR MSTR
	       :MSE MSE
	       :F* F*
	       )))

;;; Tukey's HSD post hoc, using info from http://vassun.vassar.edu/~lowry/ch14pt2.html
;;; Since the number of observations is always the same, we don't have to worry about
;;; using harminoc means, or any such gunk.  This does all the pairwise comparisons,
;;; between the means, telling us which comparisons are valid.  You give it the actual
;;; values and the within group MS Error, and it does the rest.  (Note that the number
;;; of treatments is (length data), which is called K in the HSD calculation, and 
;;; = the between df+1.) What you get back is each pairwise mean, and its significance
;;; level, as n/s * or **

(defun Anova1+TukeyHSD (data)
  (let* ((result (anova1 data))
	(k (getf result :r))
	(dfwg (getf result :SSEdf))
	(partial (sqrt (/ (getf result :mse) (length (car data))))) ; cheat -- assumes same n for all -- should use harmonic mean
	(q2 (loop for (dflimit . p2) in (let ((q (reverse (copy-list (cdr (assoc k q-table)))))) ;  UUU wants to use the next SMALLER value...
					  (or q ; Maybe not in table-- just use the highest (last) entry)
					      (reverse (copy-list (cdr (car (last q-table)))))))
		 until (<= dflimit dfwg)
		 finally (return p2)))
	(q05 (first q2))
	(q01 (second q2))
	)
    (setf (getf result :post-hocs-by-HSD-05) (TukeyHSD (getf result :meanYi-) (* q05 partial)))
    (setf (getf result :post-hocs-by-HSD-01) (TukeyHSD (getf result :meanYi-) (* q01 partial)))
    result))

(defun TukeyHSD (means q)
  (cond ((null means) ())
	(t (append (mapcar #'(lambda (e) (list (car means) e 
					       (if (> (abs (- (car means) e)) q) '+ '-))) 
			   (cdr means))
		   (TukeyHSD (cdr means) q)))))

;;; This is a convenience that lets you summarize the TukeyHSD results in a slightly more
;;; legible form.   It takes the plist that results from Anova1+TukeyHSD and gives you
;;; a printable form for the post hocs.

(defun pretty-TukeyHSD (result)
  (let* ((p05 (getf result :post-hocs-by-HSD-05))
	 (p01 (getf result :post-hocs-by-HSD-01))
	 )
    (loop for (m105 m205 +-05) in p05
	  as  (m101 m201 +-01) in p01
	  when (or (eq '+ +-05) (eq '+ +-01))
	  collect (list m105 m205 
			(if (eq '+ +-01) '** '*)))))

(setq q-table 
  '(
    (3 ; k=3
     (9     3.95   5.43)
     (10     3.88   5.27)
     (11     3.82   5.15)
     (12     3.77   5.05)
     (13     3.73   4.96)
     (14     3.70   4.89)
     (15     3.67   4.84)
     (16     3.65   4.79)
     (17     3.63   4.74)
     (18     3.61   4.70)
     (19     3.59   4.67)
     (20     3.58   4.64)
     (22     3.55   4.59)
     (24     3.53   4.55)
     (27     3.50   4.49)
     (30     3.49   4.45)
     (35     3.46   4.39)
     (40     3.44   4.37)
     (50     3.41   4.32)
     (60     3.40   4.28)
     (90     3.37   4.23)
     )
    (4 ; K=4
     (9     4.41   5.96)
     (10     4.33   5.77)
     (11     4.26   5.62)
     (12     4.20   5.50)
     (13     4.15   5.40)
     (14     4.11   5.32)
     (15     4.08   5.25)
     (16     4.05   5.19)
     (17     4.02   5.14)
     (18     4.00   5.09)
     (19     3.98   5.05)
     (20     3.96   5.02)
     (22     3.92   4.96)
     (24     3.90   4.91)
     (27     3.87   4.85)
     (30     3.85   4.80)
     (35     3.81   4.74)
     (40     3.79   4.70)
     (50     3.76   4.63)
     (60     3.74   4.59)
     (90     3.70   4.54)
     )
    (5
     (9     4.76   6.35)
     (10     4.65   6.14)
     (11     4.57   5.97)
     (12     4.51   5.84)
     (13     4.45   5.73)
     (14     4.41   5.63)
     (15     4.37   5.56)
     (16     4.33   5.49)
     (17     4.30   5.43)
     (18     4.28   5.38)
     (19     4.25   5.33)
     (20     4.23   5.29)
     (22     4.19   5.22)
     (24     4.17   5.17)
     (27     4.13   5.10)
     (30     4.10   5.05)
     (35     4.06   4.98)
     (40     4.04   4.93)
     (50     4.00   4.87)
     (60     3.98   4.82)
     (90     3.94   4.76)
     )
    ))

;;; --- Simple linear regression.

(defun regress (x y)
  (let* ((n (float (length x)))
	 (sumx (sum x))
	 (sumy (sum y))
	 (sumxy (sum (mapcar #'* x y)))
	 (sumx2 (sum (mapcar #'* x x)))
	 (m (/ (- (* n sumxy) (* sumx sumy))
	       (- (* n sumx2) (expt sumx 2))))
	 (b (+ (/ (* (- m) sumx) n)
	       (/ sumy n)))
	 (sumy2 (sum (mapcar #'* y y)))
	 (resids (mapcar #'(lambda (x y) (- y (+ (* x m) b))) x y))
	 (r (/ (- (* n sumxy) (* sumx sumy))
	       (sqrt (* 
		      (- (* n sumx2) (expt (sum x) 2))
		      (- (* n sumy2) (expt (sum y) 2))
		      ))))
	 (r2 (expt r 2))
	 )
    (list :m m :b b :resids resids :r r :r2 r2)))

;;; --- Correlation of two sequences, as in Ferguson & Takane, 1989,
;;; p. 125.  Assumes NO MISSING VALUES!

(defun correlate (x y)
  (if (not (= (length x) (length y)))
      (break "Can only correlate equal-sized sets."))
  (let* ((mx (mean x))
         (my (mean y))
	 (n (length x))
         (devx (mapcar #'(lambda (v) (- v mx)) x))
         (devy (mapcar #'(lambda (v) (- v my)) y))
	 (sumdevxy (sum (mapcar #'* devx devy)))
	 (sumsqdevx (sum (sqr devx)))
	 (sumsqdevy (sum (sqr devy)))
	 (r (/ sumdevxy (sqrt (* sumsqdevx sumsqdevy))))
	 )
    (list :r r :r2 (sqr r) :n n :p (2-tailed-correlation-significance n (abs r)))
    ))

(defun 2-tailed-correlation-significance (n r)
  ;; We use the first line for anything less than 5, and the last line for anything over 500
  ;; Otherwise, find the nearest value (maybe we should interpolate ... too much bother!)
  (let ((target-row (first *critical-values-of-r*)))
    (when (> n 5)
	  (loop for row in *critical-values-of-r*
		as (row-n . ignore) = row
		with last-row-n = 5
		with last-row = target-row
		do
		(cond ((= row-n n)
		       (setq target-row row)
		       (return t))
		      ((> row-n n)
		       (cond ((< (abs (- n row-n))
				 (abs (- n last-row-n)))
			      (setq target-row row))
			     (t (setq target-row last-row)))
		       (return t)))
		(setq last-row row)
		(setq last-row-n row-n)
		finally (progn (setq target-row row)
			       (return t))))
    (pop target-row) ; removes the N header
    (cond ((< r (car target-row)) ">0.2")
	  (t (loop for crit in (cdr target-row)
		   as p in (cdr *critical-values-of-r-two-tailed-column-interpretaion*)
		   with last-p = (car *critical-values-of-r-two-tailed-column-interpretaion*)
		   when (< r crit)
		   do (return (format nil "<~a" last-p))
		   else do (setq last-p p)
		   finally (return (format nil "<~a" p))
		   )))
    ))

;;;  Critical Values of r

;;; One tailed is half of this:
(defvar *critical-values-of-r-two-tailed-column-interpretaion* '(0.2 0.1 0.05 0.02 0.01 0.001))

(defvar *critical-values-of-r* '(
; n ... 2-tailed testing / (1-tailed testing)
;   0.2 (0.1) 0.1 (0.05) 0.05 (0.025) 0.02 (0.01) 0.01 (0.005) 0.001 (0.0005)
(5 0.687 0.805 0.878 0.934 0.959 0.991)
(6 0.608 0.729 0.811 0.882 0.917 0.974)
(7 0.551 0.669 0.754 0.833 0.875 0.951)
(8 0.507 0.621 0.707 0.789 0.834 0.925)
(9 0.472 0.582 0.666 0.750 0.798 0.898)
(10 0.443 0.549 0.632 0.715 0.765 0.872)
(11 0.419 0.521 0.602 0.685 0.735 0.847)
(12 0.398 0.497 0.576 0.658 0.708 0.823)
(13 0.380 0.476 0.553 0.634 0.684 0.801)
(14 0.365 0.458 0.532 0.612 0.661 0.780)
(15 0.351 0.441 0.514 0.592 0.641 0.760)
(16 0.338 0.426 0.497 0.574 0.623 0.742)
(17 0.327 0.412 0.482 0.558 0.606 0.725)
(18 0.317 0.400 0.468 0.543 0.590 0.708)
(19 0.308 0.389 0.456 0.529 0.575 0.693)
(20 0.299 0.378 0.444 0.516 0.561 0.679)
(21 0.291 0.369 0.433 0.503 0.549 0.665)
(22 0.284 0.360 0.423 0.492 0.537 0.652)
(23 0.277 0.352 0.413 0.482 0.526 0.640)
(24 0.271 0.344 0.404 0.472 0.515 0.629)
(25 0.265 0.337 0.396 0.462 0.505 0.618)
(26 0.260 0.330 0.388 0.453 0.496 0.607)
(27 0.255 0.323 0.381 0.445 0.487 0.597)
(28 0.250 0.317 0.374 0.437 0.479 0.588)
(29 0.245 0.311 0.367 0.430 0.471 0.579)
(30 0.241 0.306 0.361 0.423 0.463 0.570)
(40 0.207 0.264 0.312 0.367 0.403 0.501)
(50 0.184 0.235 0.279 0.328 0.361 0.451)
(60 0.168 0.214 0.254 0.300 0.330 0.414)
(80 0.145 0.185 0.220 0.260 0.286 0.361)
(100 0.129 0.165 0.197 0.232 0.256 0.324)
(120 0.118 0.151 0.179 0.212 0.234 0.297)
(140 0.109 0.140 0.166 0.196 0.217 0.275)
(160 0.102 0.130 0.155 0.184 0.203 0.258)
(180 0.096 0.123 0.146 0.173 0.192 0.243)
(200 0.091 0.117 0.139 0.164 0.182 0.231)
(300 0.074 0.095 0.113 0.134 0.149 0.189)
(400 0.064 0.082 0.098 0.116 0.129 0.164)
(500 0.057 0.074 0.088 0.104 0.115 0.147)
))

(defun even-power-of-two? (n)
  (zerop (mod (/ (log n) (log 2)) 1)))

;;; Normalize a vector by dividing it through by subtracting its min
;;; and then dividing through by its range (max-min).  If the numbers
;;; are all the same, this would screw up, so we check that first and
;;; just return a long list of 0.5 if so!

(defun normalize (v)
  (let* ((mx (mymax v))
	 (mn (mymin v))
	 (range (float (- mx mn)))
	 )
    (mapcar #'(lambda (i) (if (zerop range) 0.5
			      (/ (- i mn) range))) v)
    ))

;;; A dumb terminal way of plotting data.

(defun dumplot (v &optional show-values)
  (let* ((d1 (normalize v))
	 (d2 (mapcar #'(lambda (i) (* 50.0 i)) d1))
	 (min (mymin v))
	 (max (mymax v))
	 )
    (format t "~a~50t~a~%" min max)
    (dolist (i d2)
      (dotimes (k (round i))
        (format t " ")
	)
      (if show-values
	  (format t "* = ~a~%" (p2 (car v)))
	  (format t "*~%")
	  )
      (pop v) ; follow along for showing values
      )
    ))

;;; Cross mean takes a list of lists, as ((1 2 3) (4 3 2 1) ...) and
;;; produces a list with mean and standard error for each VERTICLE
;;; entry, so, as: ((2.5 . 1) ...) where the first pair is computed
;;; from the nth 1 of all the sublists in the input set, etc.  This is
;;; useful in some cases of data cruching.  Note that missing data is
;;; assumed to be always at the END of lists.  If it isn't, you've
;;; got to do something previously to interpolate.

(defun cross-mean (l &aux k r)
  (let* ((nmax (mymax (mapcar #'length l)))
	 (vs (make-array nmax)))
    (dolist (i l)
      (setq k 0)
      (dolist (v i)
	(push v (aref vs k))
	(incf k)))
    (dotimes (i nmax)
      (push (cons (mean (aref vs i))
		  (standard-error (aref vs i)))
	    r))
    (reverse r)))

;;; Macro to protect from division by zero.

(defmacro z/protect (expr testvar)
  `(if (zerop ,testvar) "[/0!]" ,expr))

;;; Take a set of values and produce a histogram binned into n groups, so that 
;;; you can get a report of the distribution of values.  There's a large 
;;; chance for off-by-one errores here!

(defun histovalues (v* &key (nbins 10))
  (let* ((min (min* v*))
	 (max (max* v*))
	 (inc (round (/ (abs (- min max)) nbins)))
	 (bins (loop with i = min
		     for n from 1 to nbins
		     collect (list i (incf i inc) 0)
		     ))
	 )
    (loop for v in v*
	  as bin = (loop for bin in bins
			 if (and (>= v (first bin))
				 (< v (second bin)))
			 do (incf (third bin))))
    bins))
  
;;; Simple Chi-Square

;;; This is the most basic calculation. 

(defun x2test ()
  ;; From Clarke & Cooke p. 431; should = ~7.0
  (chi-square-1 '(100 100 100 100 100 100 100 100 100 100)
		'(106  88  97 101  92 103  96 112 114  91)))

(defun chi-square-1 (expected observed)
  `(:x2 ,(loop for e in expected as o in observed
	       sum (/ (expt (- o e) 2) (float e)))
    :df ,(1- (length expected))))
  
;;; I'm not sure what the setup is supposed to be for this one, 
;;; since, like a moron I didn't give an example....

(defun chi-square-2 (table)
  (let* ((row-mars (mapcar #'(lambda (row) (apply #'+ row)) table))
	 (col-mars (loop for col from 0 to (1- (length (car table)))
			 collect (apply #'+ (loop for row in table
						  collect (nth col row)))))
	 (total (float (apply #'+ row-mars)))
	 (expectable (loop for row in table
			   as rowmar in row-mars
			   collect (loop for col in row
					 as colmar in col-mars
					 collect (cons col (/ (* rowmar colmar) total)))))
	 )
    (loop for row in expectable
	  sum (loop for entry in row
		    sum (/ (expt (- (car entry) (cdr entry)) 2) (cdr entry))))))

    
;;; This are the F score limit tables for anovas in the form: F(A,B)
;;; (From http://www.itl.nist.gov/div898/handbook/eda/section3/eda3673.htm)

;;; I *think* that A is the column df = number of columns - 1
;;;                B is the row df = number of samples - number of columns

(defun f-score>p-limit? (df1 df2 f-score limits-table)
  (let ((limit (nth df1 (assoc df2 limits-table))))
    (cond (limit (> f-score limit))
	  (t (format t "Warning! F-score>p-limit? can't find an entry for F(~a,~a)!~%" df1 df2)))))


(defvar *F0.05* '(
;   A:         1       2       3       4       5       6       7       8       9      10
;B:
(  1      161.448 199.500 215.707 224.583 230.162 233.986 236.768 238.882 240.543 241.882)
(  2       18.513  19.000  19.164  19.247  19.296  19.330  19.353  19.371  19.385  19.396)
(  3       10.128   9.552   9.277   9.117   9.013   8.941   8.887   8.845   8.812   8.786)
(  4        7.709   6.944   6.591   6.388   6.256   6.163   6.094   6.041   5.999   5.964)
(  5        6.608   5.786   5.409   5.192   5.050   4.950   4.876   4.818   4.772   4.735)
(  6        5.987   5.143   4.757   4.534   4.387   4.284   4.207   4.147   4.099   4.060)
(  7        5.591   4.737   4.347   4.120   3.972   3.866   3.787   3.726   3.677   3.637)
(  8        5.318   4.459   4.066   3.838   3.687   3.581   3.500   3.438   3.388   3.347)
(  9        5.117   4.256   3.863   3.633   3.482   3.374   3.293   3.230   3.179   3.137)
( 10        4.965   4.103   3.708   3.478   3.326   3.217   3.135   3.072   3.020   2.978)
( 11        4.844   3.982   3.587   3.357   3.204   3.095   3.012   2.948   2.896   2.854)
( 12        4.747   3.885   3.490   3.259   3.106   2.996   2.913   2.849   2.796   2.753)
( 13        4.667   3.806   3.411   3.179   3.025   2.915   2.832   2.767   2.714   2.671)
( 14        4.600   3.739   3.344   3.112   2.958   2.848   2.764   2.699   2.646   2.602)
( 15        4.543   3.682   3.287   3.056   2.901   2.790   2.707   2.641   2.588   2.544)
( 16        4.494   3.634   3.239   3.007   2.852   2.741   2.657   2.591   2.538   2.494)
( 17        4.451   3.592   3.197   2.965   2.810   2.699   2.614   2.548   2.494   2.450)
( 18        4.414   3.555   3.160   2.928   2.773   2.661   2.577   2.510   2.456   2.412)
( 19        4.381   3.522   3.127   2.895   2.740   2.628   2.544   2.477   2.423   2.378)
( 20        4.351   3.493   3.098   2.866   2.711   2.599   2.514   2.447   2.393   2.348)
( 21        4.325   3.467   3.072   2.840   2.685   2.573   2.488   2.420   2.366   2.321)
( 22        4.301   3.443   3.049   2.817   2.661   2.549   2.464   2.397   2.342   2.297)
( 23        4.279   3.422   3.028   2.796   2.640   2.528   2.442   2.375   2.320   2.275)
( 24        4.260   3.403   3.009   2.776   2.621   2.508   2.423   2.355   2.300   2.255)
( 25        4.242   3.385   2.991   2.759   2.603   2.490   2.405   2.337   2.282   2.236)
( 26        4.225   3.369   2.975   2.743   2.587   2.474   2.388   2.321   2.265   2.220)
( 27        4.210   3.354   2.960   2.728   2.572   2.459   2.373   2.305   2.250   2.204)
( 28        4.196   3.340   2.947   2.714   2.558   2.445   2.359   2.291   2.236   2.190)
( 29        4.183   3.328   2.934   2.701   2.545   2.432   2.346   2.278   2.223   2.177)
( 30        4.171   3.316   2.922   2.690   2.534   2.421   2.334   2.266   2.211   2.165)
( 31        4.160   3.305   2.911   2.679   2.523   2.409   2.323   2.255   2.199   2.153)
( 32        4.149   3.295   2.901   2.668   2.512   2.399   2.313   2.244   2.189   2.142)
( 33        4.139   3.285   2.892   2.659   2.503   2.389   2.303   2.235   2.179   2.133)
( 34        4.130   3.276   2.883   2.650   2.494   2.380   2.294   2.225   2.170   2.123)
( 35        4.121   3.267   2.874   2.641   2.485   2.372   2.285   2.217   2.161   2.114)
( 36        4.113   3.259   2.866   2.634   2.477   2.364   2.277   2.209   2.153   2.106)
( 37        4.105   3.252   2.859   2.626   2.470   2.356   2.270   2.201   2.145   2.098)
( 38        4.098   3.245   2.852   2.619   2.463   2.349   2.262   2.194   2.138   2.091)
( 39        4.091   3.238   2.845   2.612   2.456   2.342   2.255   2.187   2.131   2.084)
( 40        4.085   3.232   2.839   2.606   2.449   2.336   2.249   2.180   2.124   2.077)
( 41        4.079   3.226   2.833   2.600   2.443   2.330   2.243   2.174   2.118   2.071)
( 42        4.073   3.220   2.827   2.594   2.438   2.324   2.237   2.168   2.112   2.065)
( 43        4.067   3.214   2.822   2.589   2.432   2.318   2.232   2.163   2.106   2.059)
( 44        4.062   3.209   2.816   2.584   2.427   2.313   2.226   2.157   2.101   2.054)
( 45        4.057   3.204   2.812   2.579   2.422   2.308   2.221   2.152   2.096   2.049)
( 46        4.052   3.200   2.807   2.574   2.417   2.304   2.216   2.147   2.091   2.044)
( 47        4.047   3.195   2.802   2.570   2.413   2.299   2.212   2.143   2.086   2.039)
( 48        4.043   3.191   2.798   2.565   2.409   2.295   2.207   2.138   2.082   2.035)
( 49        4.038   3.187   2.794   2.561   2.404   2.290   2.203   2.134   2.077   2.030)
( 50        4.034   3.183   2.790   2.557   2.400   2.286   2.199   2.130   2.073   2.026)
( 51        4.030   3.179   2.786   2.553   2.397   2.283   2.195   2.126   2.069   2.022)
( 52        4.027   3.175   2.783   2.550   2.393   2.279   2.192   2.122   2.066   2.018)
( 53        4.023   3.172   2.779   2.546   2.389   2.275   2.188   2.119   2.062   2.015)
( 54        4.020   3.168   2.776   2.543   2.386   2.272   2.185   2.115   2.059   2.011)
( 55        4.016   3.165   2.773   2.540   2.383   2.269   2.181   2.112   2.055   2.008)
( 56        4.013   3.162   2.769   2.537   2.380   2.266   2.178   2.109   2.052   2.005)
( 57        4.010   3.159   2.766   2.534   2.377   2.263   2.175   2.106   2.049   2.001)
( 58        4.007   3.156   2.764   2.531   2.374   2.260   2.172   2.103   2.046   1.998)
( 59        4.004   3.153   2.761   2.528   2.371   2.257   2.169   2.100   2.043   1.995)
( 60        4.001   3.150   2.758   2.525   2.368   2.254   2.167   2.097   2.040   1.993)
( 61        3.998   3.148   2.755   2.523   2.366   2.251   2.164   2.094   2.037   1.990)
( 62        3.996   3.145   2.753   2.520   2.363   2.249   2.161   2.092   2.035   1.987)
( 63        3.993   3.143   2.751   2.518   2.361   2.246   2.159   2.089   2.032   1.985)
( 64        3.991   3.140   2.748   2.515   2.358   2.244   2.156   2.087   2.030   1.982)
( 65        3.989   3.138   2.746   2.513   2.356   2.242   2.154   2.084   2.027   1.980)
( 66        3.986   3.136   2.744   2.511   2.354   2.239   2.152   2.082   2.025   1.977)
( 67        3.984   3.134   2.742   2.509   2.352   2.237   2.150   2.080   2.023   1.975)
( 68        3.982   3.132   2.740   2.507   2.350   2.235   2.148   2.078   2.021   1.973)
( 69        3.980   3.130   2.737   2.505   2.348   2.233   2.145   2.076   2.019   1.971)
( 70        3.978   3.128   2.736   2.503   2.346   2.231   2.143   2.074   2.017   1.969)
( 71        3.976   3.126   2.734   2.501   2.344   2.229   2.142   2.072   2.015   1.967)
( 72        3.974   3.124   2.732   2.499   2.342   2.227   2.140   2.070   2.013   1.965)
( 73        3.972   3.122   2.730   2.497   2.340   2.226   2.138   2.068   2.011   1.963)
( 74        3.970   3.120   2.728   2.495   2.338   2.224   2.136   2.066   2.009   1.961)
( 75        3.968   3.119   2.727   2.494   2.337   2.222   2.134   2.064   2.007   1.959)
( 76        3.967   3.117   2.725   2.492   2.335   2.220   2.133   2.063   2.006   1.958)
( 77        3.965   3.115   2.723   2.490   2.333   2.219   2.131   2.061   2.004   1.956)
( 78        3.963   3.114   2.722   2.489   2.332   2.217   2.129   2.059   2.002   1.954)
( 79        3.962   3.112   2.720   2.487   2.330   2.216   2.128   2.058   2.001   1.953)
( 80        3.960   3.111   2.719   2.486   2.329   2.214   2.126   2.056   1.999   1.951)
( 81        3.959   3.109   2.717   2.484   2.327   2.213   2.125   2.055   1.998   1.950)
( 82        3.957   3.108   2.716   2.483   2.326   2.211   2.123   2.053   1.996   1.948)
( 83        3.956   3.107   2.715   2.482   2.324   2.210   2.122   2.052   1.995   1.947)
( 84        3.955   3.105   2.713   2.480   2.323   2.209   2.121   2.051   1.993   1.945)
( 85        3.953   3.104   2.712   2.479   2.322   2.207   2.119   2.049   1.992   1.944)
( 86        3.952   3.103   2.711   2.478   2.321   2.206   2.118   2.048   1.991   1.943)
( 87        3.951   3.101   2.709   2.476   2.319   2.205   2.117   2.047   1.989   1.941)
( 88        3.949   3.100   2.708   2.475   2.318   2.203   2.115   2.045   1.988   1.940)
( 89        3.948   3.099   2.707   2.474   2.317   2.202   2.114   2.044   1.987   1.939)
( 90        3.947   3.098   2.706   2.473   2.316   2.201   2.113   2.043   1.986   1.938)
( 91        3.946   3.097   2.705   2.472   2.315   2.200   2.112   2.042   1.984   1.936)
( 92        3.945   3.095   2.704   2.471   2.313   2.199   2.111   2.041   1.983   1.935)
( 93        3.943   3.094   2.703   2.470   2.312   2.198   2.110   2.040   1.982   1.934)
( 94        3.942   3.093   2.701   2.469   2.311   2.197   2.109   2.038   1.981   1.933)
( 95        3.941   3.092   2.700   2.467   2.310   2.196   2.108   2.037   1.980   1.932)
( 96        3.940   3.091   2.699   2.466   2.309   2.195   2.106   2.036   1.979   1.931)
( 97        3.939   3.090   2.698   2.465   2.308   2.194   2.105   2.035   1.978   1.930)
( 98        3.938   3.089   2.697   2.465   2.307   2.193   2.104   2.034   1.977   1.929)
( 99        3.937   3.088   2.696   2.464   2.306   2.192   2.103   2.033   1.976   1.928)
(100        3.936   3.087   2.696   2.463   2.305   2.191   2.103   2.032   1.975   1.927)
))

(defvar *F0.10* '(

;A:   1      2      3      4      5      6      7      8      9      10
;B:   
(  1       39.863  49.500  53.593  55.833  57.240  58.204  58.906  59.439  59.858  60.195)
(  2        8.526   9.000   9.162   9.243   9.293   9.326   9.349   9.367   9.381   9.392)
(  3        5.538   5.462   5.391   5.343   5.309   5.285   5.266   5.252   5.240   5.230)
(  4        4.545   4.325   4.191   4.107   4.051   4.010   3.979   3.955   3.936   3.920)
(  5        4.060   3.780   3.619   3.520   3.453   3.405   3.368   3.339   3.316   3.297)
(  6        3.776   3.463   3.289   3.181   3.108   3.055   3.014   2.983   2.958   2.937)
(  7        3.589   3.257   3.074   2.961   2.883   2.827   2.785   2.752   2.725   2.703)
(  8        3.458   3.113   2.924   2.806   2.726   2.668   2.624   2.589   2.561   2.538)
(  9        3.360   3.006   2.813   2.693   2.611   2.551   2.505   2.469   2.440   2.416)
( 10        3.285   2.924   2.728   2.605   2.522   2.461   2.414   2.377   2.347   2.323)
( 11        3.225   2.860   2.660   2.536   2.451   2.389   2.342   2.304   2.274   2.248)
( 12        3.177   2.807   2.606   2.480   2.394   2.331   2.283   2.245   2.214   2.188)
( 13        3.136   2.763   2.560   2.434   2.347   2.283   2.234   2.195   2.164   2.138)
( 14        3.102   2.726   2.522   2.395   2.307   2.243   2.193   2.154   2.122   2.095)
( 15        3.073   2.695   2.490   2.361   2.273   2.208   2.158   2.119   2.086   2.059)
( 16        3.048   2.668   2.462   2.333   2.244   2.178   2.128   2.088   2.055   2.028)
( 17        3.026   2.645   2.437   2.308   2.218   2.152   2.102   2.061   2.028   2.001)
( 18        3.007   2.624   2.416   2.286   2.196   2.130   2.079   2.038   2.005   1.977)
( 19        2.990   2.606   2.397   2.266   2.176   2.109   2.058   2.017   1.984   1.956)
( 20        2.975   2.589   2.380   2.249   2.158   2.091   2.040   1.999   1.965   1.937)
( 21        2.961   2.575   2.365   2.233   2.142   2.075   2.023   1.982   1.948   1.920)
( 22        2.949   2.561   2.351   2.219   2.128   2.060   2.008   1.967   1.933   1.904)
( 23        2.937   2.549   2.339   2.207   2.115   2.047   1.995   1.953   1.919   1.890)
( 24        2.927   2.538   2.327   2.195   2.103   2.035   1.983   1.941   1.906   1.877)
( 25        2.918   2.528   2.317   2.184   2.092   2.024   1.971   1.929   1.895   1.866)
( 26        2.909   2.519   2.307   2.174   2.082   2.014   1.961   1.919   1.884   1.855)
( 27        2.901   2.511   2.299   2.165   2.073   2.005   1.952   1.909   1.874   1.845)
( 28        2.894   2.503   2.291   2.157   2.064   1.996   1.943   1.900   1.865   1.836)
( 29        2.887   2.495   2.283   2.149   2.057   1.988   1.935   1.892   1.857   1.827)
( 30        2.881   2.489   2.276   2.142   2.049   1.980   1.927   1.884   1.849   1.819)
( 31        2.875   2.482   2.270   2.136   2.042   1.973   1.920   1.877   1.842   1.812)
( 32        2.869   2.477   2.263   2.129   2.036   1.967   1.913   1.870   1.835   1.805)
( 33        2.864   2.471   2.258   2.123   2.030   1.961   1.907   1.864   1.828   1.799)
( 34        2.859   2.466   2.252   2.118   2.024   1.955   1.901   1.858   1.822   1.793)
( 35        2.855   2.461   2.247   2.113   2.019   1.950   1.896   1.852   1.817   1.787)
( 36        2.850   2.456   2.243   2.108   2.014   1.945   1.891   1.847   1.811   1.781)
( 37        2.846   2.452   2.238   2.103   2.009   1.940   1.886   1.842   1.806   1.776)
( 38        2.842   2.448   2.234   2.099   2.005   1.935   1.881   1.838   1.802   1.772)
( 39        2.839   2.444   2.230   2.095   2.001   1.931   1.877   1.833   1.797   1.767)
( 40        2.835   2.440   2.226   2.091   1.997   1.927   1.873   1.829   1.793   1.763)
( 41        2.832   2.437   2.222   2.087   1.993   1.923   1.869   1.825   1.789   1.759)
( 42        2.829   2.434   2.219   2.084   1.989   1.919   1.865   1.821   1.785   1.755)
( 43        2.826   2.430   2.216   2.080   1.986   1.916   1.861   1.817   1.781   1.751)
( 44        2.823   2.427   2.213   2.077   1.983   1.913   1.858   1.814   1.778   1.747)
( 45        2.820   2.425   2.210   2.074   1.980   1.909   1.855   1.811   1.774   1.744)
( 46        2.818   2.422   2.207   2.071   1.977   1.906   1.852   1.808   1.771   1.741)
( 47        2.815   2.419   2.204   2.068   1.974   1.903   1.849   1.805   1.768   1.738)
( 48        2.813   2.417   2.202   2.066   1.971   1.901   1.846   1.802   1.765   1.735)
( 49        2.811   2.414   2.199   2.063   1.968   1.898   1.843   1.799   1.763   1.732)
( 50        2.809   2.412   2.197   2.061   1.966   1.895   1.840   1.796   1.760   1.729)
( 51        2.807   2.410   2.194   2.058   1.964   1.893   1.838   1.794   1.757   1.727)
( 52        2.805   2.408   2.192   2.056   1.961   1.891   1.836   1.791   1.755   1.724)
( 53        2.803   2.406   2.190   2.054   1.959   1.888   1.833   1.789   1.752   1.722)
( 54        2.801   2.404   2.188   2.052   1.957   1.886   1.831   1.787   1.750   1.719)
( 55        2.799   2.402   2.186   2.050   1.955   1.884   1.829   1.785   1.748   1.717)
( 56        2.797   2.400   2.184   2.048   1.953   1.882   1.827   1.782   1.746   1.715)
( 57        2.796   2.398   2.182   2.046   1.951   1.880   1.825   1.780   1.744   1.713)
( 58        2.794   2.396   2.181   2.044   1.949   1.878   1.823   1.779   1.742   1.711)
( 59        2.793   2.395   2.179   2.043   1.947   1.876   1.821   1.777   1.740   1.709)
( 60        2.791   2.393   2.177   2.041   1.946   1.875   1.819   1.775   1.738   1.707)
( 61        2.790   2.392   2.176   2.039   1.944   1.873   1.818   1.773   1.736   1.705)
( 62        2.788   2.390   2.174   2.038   1.942   1.871   1.816   1.771   1.735   1.703)
( 63        2.787   2.389   2.173   2.036   1.941   1.870   1.814   1.770   1.733   1.702)
( 64        2.786   2.387   2.171   2.035   1.939   1.868   1.813   1.768   1.731   1.700)
( 65        2.784   2.386   2.170   2.033   1.938   1.867   1.811   1.767   1.730   1.699)
( 66        2.783   2.385   2.169   2.032   1.937   1.865   1.810   1.765   1.728   1.697)
( 67        2.782   2.384   2.167   2.031   1.935   1.864   1.808   1.764   1.727   1.696)
( 68        2.781   2.382   2.166   2.029   1.934   1.863   1.807   1.762   1.725   1.694)
( 69        2.780   2.381   2.165   2.028   1.933   1.861   1.806   1.761   1.724   1.693)
( 70        2.779   2.380   2.164   2.027   1.931   1.860   1.804   1.760   1.723   1.691)
( 71        2.778   2.379   2.163   2.026   1.930   1.859   1.803   1.758   1.721   1.690)
( 72        2.777   2.378   2.161   2.025   1.929   1.858   1.802   1.757   1.720   1.689)
( 73        2.776   2.377   2.160   2.024   1.928   1.856   1.801   1.756   1.719   1.687)
( 74        2.775   2.376   2.159   2.022   1.927   1.855   1.800   1.755   1.718   1.686)
( 75        2.774   2.375   2.158   2.021   1.926   1.854   1.798   1.754   1.716   1.685)
( 76        2.773   2.374   2.157   2.020   1.925   1.853   1.797   1.752   1.715   1.684)
( 77        2.772   2.373   2.156   2.019   1.924   1.852   1.796   1.751   1.714   1.683)
( 78        2.771   2.372   2.155   2.018   1.923   1.851   1.795   1.750   1.713   1.682)
( 79        2.770   2.371   2.154   2.017   1.922   1.850   1.794   1.749   1.712   1.681)
( 80        2.769   2.370   2.154   2.016   1.921   1.849   1.793   1.748   1.711   1.680)
( 81        2.769   2.369   2.153   2.016   1.920   1.848   1.792   1.747   1.710   1.679)
( 82        2.768   2.368   2.152   2.015   1.919   1.847   1.791   1.746   1.709   1.678)
( 83        2.767   2.368   2.151   2.014   1.918   1.846   1.790   1.745   1.708   1.677)
( 84        2.766   2.367   2.150   2.013   1.917   1.845   1.790   1.744   1.707   1.676)
( 85        2.765   2.366   2.149   2.012   1.916   1.845   1.789   1.744   1.706   1.675)
( 86        2.765   2.365   2.149   2.011   1.915   1.844   1.788   1.743   1.705   1.674)
( 87        2.764   2.365   2.148   2.011   1.915   1.843   1.787   1.742   1.705   1.673)
( 88        2.763   2.364   2.147   2.010   1.914   1.842   1.786   1.741   1.704   1.672)
( 89        2.763   2.363   2.146   2.009   1.913   1.841   1.785   1.740   1.703   1.671)
( 90        2.762   2.363   2.146   2.008   1.912   1.841   1.785   1.739   1.702   1.670)
( 91        2.761   2.362   2.145   2.008   1.912   1.840   1.784   1.739   1.701   1.670)
( 92        2.761   2.361   2.144   2.007   1.911   1.839   1.783   1.738   1.701   1.669)
( 93        2.760   2.361   2.144   2.006   1.910   1.838   1.782   1.737   1.700   1.668)
( 94        2.760   2.360   2.143   2.006   1.910   1.838   1.782   1.736   1.699   1.667)
( 95        2.759   2.359   2.142   2.005   1.909   1.837   1.781   1.736   1.698   1.667)
( 96        2.759   2.359   2.142   2.004   1.908   1.836   1.780   1.735   1.698   1.666)
( 97        2.758   2.358   2.141   2.004   1.908   1.836   1.780   1.734   1.697   1.665)
( 98        2.757   2.358   2.141   2.003   1.907   1.835   1.779   1.734   1.696   1.665)
( 99        2.757   2.357   2.140   2.003   1.906   1.835   1.778   1.733   1.696   1.664)
(100        2.756   2.356   2.139   2.002   1.906   1.834   1.778   1.732   1.695   1.663)
))

;;; Nonparametric one-sample (signed) rank test (Wilcoxon):

#|

 From :http://www.graphpad.com/instatman/HowtheWilcoxonranksumtestworks.htm

InStat follows these steps: 

1. Calculate how far each value is from the hypothetical value.

2. Ignore values that exactly equal the hypothetical value. Call the
number of remaining values N.

3. Rank these distances, paying no attention to whether the values are
higher or lower than the hypothetical value.

4. For each value that is lower than the hypothetical value, multiply
the rank by negative 1.

5. Sum the positive ranks. InStat reports this value.

6. Sum the negative ranks. InStat also reports this value.

7. Add the two sums together. This is the sum of signed ranks, which
InStat reports as W.

If the data really were sampled from a population with the
hypothetical mean, you'd expect W to be near zero. If W (the sum of
signed ranks) is far from zero, the P value will be small. The P value
answers this question: Assume that you randomly sample N values from a
population with the hypothetical median. What is the chance that W
will be as far from zero (or further) as you observed?
|#

(defun wilcoxon-1 (initial-values target)
  (let ((n 0) 
	r/v/d* 
	(sum+ranks 0)
	(sum-ranks 0))
    ;; 1. Calculate how far each value is from the hypothetical value.
    ;; 2. Ignore values that exactly equal the hypothetical value. Call the
    ;;   number of remaining values N.
    (dolist (v initial-values)
      (when (not (= target v))
	(push (list 'rank v (abs (- target v))) r/v/d*)
	(incf n)))
    ;; 3. Rank these distances, paying no attention to whether the values are
    ;;    higher or lower than the hypothetical value.
    (setq r/v/d* (sort r/v/d* #'(lambda (a b) (< (third a) (third b)))))
    ;; Ranking doesn't deal with ties!!!
    (loop for entry in r/v/d*
	  as rank from 1 by 1
	  do (setf (car entry) rank))
    ;; 4. For each value that is lower than the hypothetical value, multiply
    ;;    the rank by negative 1.
    (loop for entry in r/v/d*
	  as (rank value distance) = entry
	  when (< value target)
	  do (setf (car entry) (- rank)))
    ;; 5. Sum the positive ranks. InStat reports this value.
    ;; 6. Sum the negative ranks. InStat also reports this value.
    ;; 7. Add the two sums together. This is the sum of signed ranks, which
    ;;    InStat reports as W.
    (print r/v/d*)
    (loop for (rank value distance) in r/v/d*
	  do (cond ((> rank 0) (incf sum+ranks rank))
		   ((< rank 0) (incf sum-ranks rank))
		   )
	  finally (return (values sum+ranks sum-ranks (+ sum+ranks sum-ranks))))
    ))

#|

;;; Rank a list properly, including dealing with ties.  This takes a 
;;; value list (30 20 30 40 50 60 ...) and returns the list sorted and 
;;; ranked in pairs, as (for the above): ((20 . 1) (30 . 2.5) (30 . 2.5) (40 . 4) (50 . 5)...)
;;; Notice that the two 30's split the rank position for 2 and 3 equally; Any number 
;;; of ties will split the positional difference, so, three ties will end up on an integer.
;;; This is done in two passes.  First we sort and assign sequential integers.  Then
;;; we re-assign ties.

(defun tr ()
  (rank-values (copy-list '(80 40 50 20 60 10 30 70))))

(defun rank-values (v*)
  (let* ((v/r* (loop for v in v* collect (cons v 'rank)))
	 (v/r* (sort v/r* #'(lambda (a b) (< (car a) (car b)))))
	 (ignore (loop for v/r in v/r* as rank from 1 by 1 
		       do (setf (cdr v/r) rank)))
	 )

    (loop for v/r+ on v/r*
	  as major-rank = from 1 by 1
	  with minor-rank = 1
	  as count = (loop with target = (caar v/r+)
			   as (next . minor-rank) in (cdr v/r+)
			   as k from 1 by 1
			   until (not (= next target))
			   finally (return k))
	  do 
	  (print (list (car v/r+) minor-rank count))
	  (when (not (= 1 count))
		(let ((new-rank (+ minor-rank (/ count 2.0))))
		  (loop with target = (caar v/r+)
			as entry in (cdr v/r+)
			as (next . minor-rank) = entry
			until (not (= next target))
			do (setf (cdr entry) new-rank))))
	  (setf minor-rank (+ rank count))
	  )
    v/r*))
|#

