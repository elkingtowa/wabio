(defun whtest ()
  (let* ((test-in-words
	 '(("cat" . "dog")
	   ("elephant" . "mouse")
	   ("dinosaur" . "trex")
	   ("horse" . "cow")
	   ("snake" . "apple")
	   ("robin" . "batman")
	   ))
	 (testwh (initwh test-in-words))
	 (test-out-words
	  '("car" "house" "scan" "ribbon"))
	 )
    (loop for word in test-out-words
	  do (print (list word (sort (word-homology word testwh)
				     #'(lambda (a b) (> (car a) (car b))))))
	  )
    ))

;;; This passes around whinfo structs that contain the relevant tables so that 
;;; you can have multiple word homology datasets.

;;; lsht is the Letter sequence Hash Table
;;; lsk is the Letter sequence counter, incremented for each unique pair of letters
;;; c is the list of (word . result) where what is returned when a word is found, is result 
;;; tl is the compiled target list

(defstruct whinfo c lsht lsk tl)

;;; This function initilizes the word-homology database for a given whinfo struct

(defun initwh (c)
  (let ((whinfo (make-whinfo :c c)))
    (setf (whinfo-lsht whinfo) (make-hash-table :test #'equal))
    (setf (whinfo-lsk whinfo) 0)
    (setf (whinfo-tl whinfo) (compile-target-list (whinfo-c whinfo) whinfo))
    whinfo))

;;; This converts a word to a series of numbers where each unique pair of letters 
;;; is a new number and adds it to the existing word-homology database.  The resulting
;;; list of numbers is sorted so that fast-n-member can cut off search early.  Note
;;; that the pairs are order-independent, so that AB is the same as BA although they
;;; are only stored in one order in the table (which is the first order that they
;;; happen to be seen in).

(defun compile-word (word whinfo)
  (setq word (string-downcase word))
  (sort 
   (loop for p from 0 to (- (length word) 2)
	 as l1 = (aref word p)
	 as l2 = (aref word (1+ p))
	 as l1.l2 = (cons l1 l2)
	 as l2.l1 = (cons l2 l1)
	 as k = (or (gethash l1.l2 (whinfo-lsht whinfo))
		    (gethash l2.l1 (whinfo-lsht whinfo)))
	 collect (or k
		     (let ((k (incf (whinfo-lsk whinfo))))
		       (setf (gethash l1.l2 (whinfo-lsht whinfo)) k)
		       k)))
   #'<))

;;; This function determines how similar the two input lists are, where the lists are the numbers generated
;;; by compile-word

(defun score-homology (w1 w2)
  (let* ((o1 (loop for l1 in w1 if (fast-n-member l1 w2) sum 1))
	 (o2 (loop for l2 in w2 if (fast-n-member l2 w1) sum 1))
	 (l1 (length w1))
	 (l2 (length w2))
	 )
    (/ (+ (/ o1 l1) (/ o2 l2)) 2.0)))

;;; This does member for numbers, but assumes that the target list is sorted
;;; so that if it hits a number greater than the one we're at, it can cut the
;;; search short.

(defun fast-n-member (n n*)
  (loop for n2 in n*
	do (cond ((= n n2) (return t))
		 ((> n2 n) (return nil)))))

;;; This compiles a list of words with compile-word

(defun compile-target-list (names whinfo)
  (loop for name in names
	collect (cons (compile-word (car name) whinfo) name)))

;;; Get the top N scoring of a set of words w/o having to sort, which is too slow.

(defun word-homology (w whinfo &optional (n 3))
  (let* ((wc (compile-word w whinfo))
	 (topn (loop for i from 1 to n collect (list 0 'foo 'bar)))
	 )
    (when wc ; Compile-word gives nil if there's no entries that match any letter combination
	  (if (not (whinfo-tl whinfo)) (import-kb))
	  ;; If the new word is more than the lowest of the ones in the set, replace that one with the new one.
	  (loop for (tc . tw) in (whinfo-tl whinfo)
		as score = (score-homology wc tc)
		as lowest = (loop with least-entry = (car topn)
				  with least-score = (car least-entry)
				  for test-entry in (cdr topn)
				  as ts = (car test-entry)
				  if (< ts least-score)
				  do 
				  (setq least-score ts least-entry test-entry)
				  finally (return least-entry))
		if (> score (car lowest))
		do 
		(setf (first lowest) score)
		(setf (second lowest) tw)
		(setf (third lowest) w)
		finally (return (remove (list 0 'foo 'bar) topn :test #'equal))
		))
    ))

