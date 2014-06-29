;;; Common Lisp Utilities Copyright (c) Jeff Shrager 1999-2006 
;;; Contact: jshrager@stanford.edu

;;; ===================================================================
;;; Various levels of reporting.  Change this to LOUD LOW or SILENT
;;; SCREAM reports get through LOW, but WHISPER reports don't.  SILENT
;;; cuts everything off.

(defvar *volume* 'loud)
(defvar *silence-counter* 0)

(defmacro scream (format &rest args)
  `(cond ((eq *volume* 'silent) (check-silence))
	 (t (format t ,format ,@args))))
(defmacro whisper (format &rest args)
  `(progn (when (not (member *volume* '(low silent)))
		(format t ,format ,@args))
	  (when (eq *volume* 'silent) (check-silence))))
  (defmacro volume (?)
  `(case ',?
	((loud high) (setq *volume* 'loud))
	((low whisper) (setq *volume* 'low))
	((off down silent) (setq *volume* 'silent))
	(t "Error: Volume must be one of LOUD/HIGH, LOW/WHISPER, or OFF/DOWN/SILENT")))
(defun check-silence ()
  (when (zerop (mod (incf *silence-counter*) 1000))
	(format t "[Warning, you have the cone of silence down!  Use: (volume low) or (volume loud) to turn it up.]~%")))

;;; Used to display the values of things for debugging
;;; Used to be called (Display ...)

(defmacro display (&rest ?)
  `(*foo* ,@?))

(defmacro *FOO* (&rest names)
  `(progn 
     (setq *FOO*
	   (list ,@(loop for name in names
		   collect `(cons ',name ,name))))
     (loop for (name . value) in *FOO*
	   do (format t "~a = ~s~%" name value))
     (cdar (last *FOO*))))

(defun keywordize (symbol &optional (case :upper))
  "Return a symbol in the KEYWORD package with the same name as SYMBOL (modified appropriately for CASE).  SYMBOL may also be a string."
  (let ((p (find-package :keyword)))
    (flet ((doit 
	    (string)
	    (intern
	     (ecase case
		    (:upper (string-upcase string))
		    (:lower (string-downcase string))
		    ((:preserve nil) string))
	     p)))
	  (cond
	   ((stringp symbol) (doit symbol))
	   ((symbolp symbol) (doit (symbol-name symbol)))
	   (t (error "Invalid argument to KEYWORDIZE: ~A" symbol))
	   ))))

;;; Displays the first n entries in a hash table:

(defun dht (table &optional (n 10))
  (maphash #'(lambda (key value)
	       (when (zerop (decf n)) (return-from dht))
	       (format t "~s: ~s~%" key value)	       
	       )
	   table))

;;; --- Some input utilities.

(defun pause ()
  (if (not (y-or-n-p "Continue?"))
      (break)))

(defun random-one-of (l)
  (nth (random (length l)) l))

(defun apply-across-dirtree (path fn &key per-dir-initfn)
  "Beginning at a given directory, run through the entire directory
   tree and do something specific. A context is passed to fn of the
   name of the current complete path. If the optional per-dir-initfn is
   given, it is applied to each newly found directory, with the pathname,
   before the directory is processed."
  (when per-dir-initfn (funcall per-dir-initfn path))
  (loop for file in (directory (format nil "~a/*.*" path))
	do 
	(cond ((probe-directory file)
	       (apply-across-dirtree (namestring file) 
				     fn
				     :per-dir-initfn per-dir-initfn))
	      (t (funcall fn (namestring file))))))

;;; ===================================================================
;;; String functions.

(defun parse-string (string &key (break-char #\space) (make-substitutions nil))
  "break the string up at spaces into little
  strings. These must be (read-from-string'ed...) outside if you
  want numbers or atoms.  We replace pounds and colons with * and -
  respectively so that if you do want to read-from-string the
  result, you aren't screwed.  In you'd rather not break things up
  on space, and know that it's all lisp-readable, you might rather
  use multi-rfs (multi-read-from-string).  Someday I ought to do
  this with readtables...."
  (when make-substitutions (fast-substitute string))
  (prog (pos item results break-string)
        (setq break-string (format nil "~a" break-char))
    loop
        (setq string (string-left-trim break-string string))
	(setq item (parse-string-from-string string break-char))
	(if (null item) (return (reverse results)))
	(push item results)
	(setq pos (position break-char string))
	(if (null pos) (return (reverse results)))
	(setq string (subseq string pos))
	(go loop)))

;;; This version takes a string of break chars.

(defun parse-string2 (string &key (make-substitutions nil) (break-string " "))
  (when make-substitutions (fast-substitute string))
  (prog (item results)
    loop
        (setq string (string-left-trim break-string string))
	(setq item (parse-string-from-string2 string break-string))
	(if (null item) (return (reverse results)))
	(push item results)
	(setq string (subseq string (length item)))
	(go loop)))

(defun parse-string-from-string2 (instring break-string)
  (if (zerop (length instring))
      ()
    (prog (outstring)
	  (setq outstring "")
      loop
          (if (or (zerop (length instring))
		  (position (aref instring 0) break-string))
	      (return outstring))
	  (setq outstring (format nil "~a~a" outstring (aref instring 0)))
	  (setq instring (subseq instring 1))
	  (go loop)
	  )
    ))

;;; WARNING: THIS DESTROYS THE STRING (which is usually an okay thing to do)
;;; (Which of these is faster? -- looks like they take *exactly* the same amount of time!)

(defun fast-substitute (string)
  (loop for c across string
	as p from 0 by 1
	when (or (char-equal c #\.)
		 (char-equal c #\,)
		 (char-equal c #\()
		 (char-equal c #\))
		 (char-equal c #\')
		 (char-equal c #\#)
		 (char-equal c #\tab)
		 (char-equal c #\:)
		 (char-equal c #\!)
		 (char-equal c #\")
		 )
	do (setf (aref string p) #\space))
  string)

(defun fast-substitute (string)
  (dotimes (p (length string))
    (let ((c (aref string p)))
      (if (char-equal c #\tab)
	  (setf (aref string p) #\space)
	(if (char-equal c #\:)
	    (setf (aref string p) #\space)
	  (if (char-equal c #\#)
	      (setf (aref string p) #\space)
	    (if (char-equal c #\()
		(setf (aref string p) #\space)
	      (if (char-equal c #\))
		  (setf (aref string p) #\space)
		(if (char-equal c #\")
		    (setf (aref string p) #\space)
		  (if (char-equal c #\')
		      (setf (aref string p) #\space)
		    (if (char-equal c #\!)
			(setf (aref string p) #\space)
		      ))))))))))
  string)

;;; Stops on space, and returns nil when there's nothing there to begin with.

(defun parse-string-from-string (instring break-char)
  (if (zerop (length instring))
      ()
    (prog (outstring)
	  (setq outstring "")
      loop
          (if (or (zerop (length instring))
		  (eq break-char (aref instring 0)))
	      (return outstring))
	  (setq outstring (format nil "~a~a" outstring (aref instring 0)))
	  (setq instring (subseq instring 1))
	  (go loop)
	  )
    ))

;;; Mike Traver's highly efficient, although somewhat less featureful
;;; version, allows string sharing thus saving conses:

;;; Split string into substrings delimited by delimiter.  Unless copy
;;; is t, return a list of substrings that share structure with the
;;; original string.  If num-values is provided, it is an integer
;;; representing the number of items to return.  Results will be
;;; padded or truncated.

(defun string-split (string &key (delimiter #\space) (copy t) num-values)
  (let ((substrings '())
        (length (length string))
        (last 0))
    (flet ((add-substring (i)
             (push (if copy
                       (subseq string last i)
                     (make-array (- i last)
                                 :element-type 'character
                                 :displaced-to string
                                 :displaced-index-offset last))
              substrings)))
      (dotimes (i length)
        (when (eq (char string i) delimiter)
          (add-substring i)
          (setq last (1+ i))))
      (add-substring length)
      (if num-values
          (pad-or-truncate num-values (nreverse substrings))
        (nreverse substrings)))))

;;; Do read-from-strings until the string is empty.  This is often used 
;;; in place of:
;;;  (mapcar #'read-from-string (parse-string ...))
;;; when you know that there's no crap in the entry.

(defun multi-read-from-string (string &aux result (start 0))
  (loop 
    (multiple-value-bind (item nextloc) 
              (read-from-string string nil '*eof* :start start)
      (if (eq '*eof* item) 
	  (return (reverse result))
	  (progn (push item result)
		 (setq start nextloc))
	  ))))

;;; ===================================================================
;;; Other stuff.

;;; Note that these prompting utils REQUIRE the default argument to be
;;; be given, which is a tad weird but it's that way so that I can use
;;; the rest for the format arg list because you can't do both &key
;;; and &rest.  I guess that I ought to use a macro here.

(defun prompt-for-number (form default &rest args &aux n)
  (setq form (format nil "~%~a [~a]: " form (if default default "no default")))
  (loop 
    (apply #'format (cons t (cons form args)))
    (let ((ans (read-line t nil nil)))
      (setq n (if (zerop (length ans))
		  () 
		  (read-from-string ans)))
	    )
    (if (null n)
        (if default
  	    (return default)
	    (format t "~%Sorry.  There's no default.~%"))
        (if (numberp n) 
	    (return n)
	    (format t "~%This has to be a number. Please try again.~%")
	    ))
    )
  )

;;; Note that prompt-for-string doesn't have any error checking, and
;;; will return nil if a null string is entered.  This might be a
;;; problem for callers.

(defun prompt-for-string (form &rest args)
  (apply #'format (append (list t) (list form) args))
  (let ((ans (read-line)))
    (if (zerop (length ans)) () ans)))

;;; This version is a little weird; it is just a short cut to use the
;;; default value, which just saves the caller from having to test and
;;; skip the prompt.

(defun get-string-w-default (prompt default &optional use-default)
  (let ((ans (if use-default 
		 (progn (format t "Using ~a~%" default) default)
	         (prompt-for-string
		   (format nil "~a (~a): " prompt default)))))
    (if ans ans default)))

;;;

(defun file-exists? (name)
  (let ((handle (open name :direction :input :if-does-not-exist nil)))
    (if handle
	(progn (close handle) t)
      ())))

;;; --- Some extensions.

;;; Same as system but does error checking.

(defun system! (cmdform &rest args)
  (let*( (cmd (apply #'format (append (list nil cmdform ) args)))
	 (r (system cmd)))
    (if (not (zerop r))
	(progn (format t "System call:[~a] failed, error ~a!~%" cmd r)
	       (break "Stopping for assessment!"))
      )))

;;; Order-irrelevant at the top level.

(defun set-equal (a b)
  (cond ((null a) (null b))
        ((member (car a) b :test #'equal)
	 (set-equal (cdr a) (remove (car a) b :test #'equal)))
	(t ())))

;;; Read a file into a list of lists of its lines, using read-from-string
;;; on each line.  The protected? arg only does the processing if the 
;;; first char of the line is a number.

(defun read-table-file (fn &optional protected? &aux l r)
  (with-open-file (f fn :direction :input)
     (loop (setq l (read-line f nil nil))
	   (if (null l) (return nil))
	   (if (or (not protected?)
		   (and (> (length l) 0)
			(let ((tc (char-code (aref l 0))))
			  (and (>= tc 48) (<= tc 57)))))
	       (push (multi-read-from-string l) r))
	       )
     (reverse r)))

;;; Do read-from-strings until the string is empty.  May be used
;;; in place of:  (mapcar #'read-from-string (parse-string ...))
;;; when you know that there's no crap in the entry, else use
;;; the safer: read-raw-table-file.

(defun multi-read-from-string (string &aux result (start 0))
  (loop 
    (multiple-value-bind (item nextloc) 
              (read-from-string string nil '*eof* :start start)
      (if (eq '*eof* item) 
	  (return (reverse result))
	  (progn (push item result)
		 (setq start nextloc))
	  ))))

;;; Safe version of read-table-file for tab delimited files.
;;; (Sort of slow since it uses parse-string; could be improved, but it's safer
;;;  since it reads everything as a string.)  Note: Returns everything as strings!

(defun read-raw-table-file (fn &key (break-char #\tab))
  (with-open-file (f fn :direction :input)
     (loop for line = (read-line f nil nil)
	   as k from 1 by 1
	   until (null line)
	   do (when (zerop (mod k 100)) (print k))
	   collect (parse-string line :break-char break-char 
				 :make-substitutions nil))))

;;; We can write either a column or matrix. 

(defun write-table-file (fn table &optional rounding-fn)
  (with-open-file (f fn :direction :output :if-exists :supersede)
    (dolist (line table)
      (if (numberp line) (setq line (list line)))
      (dolist (v line)
        (format f "~a " (if rounding-fn (funcall rounding-fn v) v)))
      (format f "~%")
      )
    ))

;;; Get the first n items of a list, somewhat efficiently.

(defun first-n (n l &aux r)
  (setq r (list '!))
  (do ((s l (cdr s))
       (v r (cdr v))
       (k n (1- k)))
      ((zerop k) (cdr r))
      (rplacd v (list (car s)))
      ))

;;; This return the LAST n items is a list; Note that it does NOT create
;;; a new list!

(defun last-n (n l)
  (nthcdr (- (length l) n) l))

;;; Form all combinations of items in a list, including the list itself.

(defun all-sublists (l)
  "From '(a s d) produce: ((D) NIL (S D) (S) (A D) (A) (A S D) (A S))"
  (cond ((null l) ())
	((null (cdr l)) (list l nil))
        (t (append (all-sublists (cdr l))
		   (mapcar #'(lambda (i) (cons (car l) i))
			   (all-sublists (cdr l)))))
	))

(defun all-ordered-sublists (l)
  "From '(a s d f) produce: ((A) (S) (D) (F) (A S) (S D) (D F) (A S D) (S D F) (A S D F))"
  (let ((len (length l)))
    (loop for curlen from 1 to len
	  append (loop for w+ on l
		       as k from 1 to (1+ (- len curlen))
		       collect (loop for i from 1 to curlen
				     as w in w+
				     collect w)))))

;;; ===================================================================
;;; --- Time/date functions.

(defun get-time (&optional (utime (get-universal-time)))
  (multiple-value-bind
    (second minute hour date month year dow dstp tz)
    (decode-universal-time (or utime ))
    (list :second second
	  :minute minute
	  :hour (+ 2 hour) ;; I don't get why this needs to be here! DST or something?
	  :date date
	  :month month
	  :year year
	  :dow dow
	  :dstp dstp
	  :tz tz)))

(defun time-as-int (&optional (utime (get-universal-time)))
  (let ((time (get-time utime)))
  (+ (* (getf time :hour) 3600)
     (getf time :second)
     (* 60 (getf time :minute)))))

(defun pretty-time (&optional (utime (get-universal-time)))
  (let* ((time-as-int (time-as-int utime))
	 (r (mod time-as-int 3600))
	 (h (truncate (/ time-as-int 3600)))
	 (s (mod r 60))
	 (m (truncate (/ r 60))))
    (format nil "~2,'0d:~2,'0d:~2,'0d" h m s)))

(defvar *months* '("Jan" "Feb" "Mar" "Apr" "May" "Jun" "Jul" "Aug" "Sep" "Oct" "Nov" "Dec"))

(defun pretty-date (&optional (utime (get-universal-time)))
  (let* ((time (get-time utime))
	 (date (getf time :date))
	 (month (elt *months* (1- (getf time :month))))
	 (year (getf time :year))
	 )
    (format nil "~a ~a ~a" date month year)))
    
(defun jeff-date () (date-as-int))

(defun date-as-int ()
  (let ((time (get-time)))
    (read-from-string
     (format nil "~4d~2,'0d~2,'0d" 
	         (getf time :year)
		 (getf time :month)
		 (getf time :date))
     )))

;;; A Jeff date can be 980731 or 19980731 (although the latter is more standard)
;;; This standardizes them.

(defun std-jeff-date (date)
  (let* ((date (format nil "~a" date)) ; ensure that it's a string
	 (len (length date)))
    (cond ((= len 6)
	   (read-from-string (format nil "19~a" date)))
	  ((= len 8)
	   (read-from-string date))
	  (t (break "Tried to standardize this jeff-date: ~a" date))
	  )))

;;; To compare Jeff Dates you need to use this function to standardize
;;; and julianize them, then you can use them as real valued days from
;;; 1/1/1900 This is only approximate since it assumes 30 days/month
;;; and 365 days/yr.

(defun norm-jeff-date (jd)
  (let* ((sjd (std-jeff-date jd))
	 (jdyear (truncate (/ sjd 10000)))
	 (jdmo (truncate (/ (- sjd (* jdyear 10000)) 100)))
	 (jdday (- sjd (+ (* jdyear 10000) (* jdmo 100))))
	 )
    (+ (* 365 (- jdyear 1900))
       (* 30 jdmo)
       jdday)))

(defvar *jdmodays* '(31 27 31 30 31 30 31 31 30 31 30 31))

(defun explode-jeff-date (jd)
  (let* ((sjd (std-jeff-date jd))
	 (jdyear (truncate (/ sjd 10000)))
	 (jdmo (truncate (/ (- sjd (* jdyear 10000)) 100)))
	 (jdday (- sjd (+ (* jdyear 10000) (* jdmo 100))))
	 )
    (list jdyear jdmo jdday)))

(defun jeff-date+1 (jd)
  (let* ((sjd (std-jeff-date jd))
	 (jdyear (truncate (/ sjd 10000)))
	 (jdmo (truncate (/ (- sjd (* jdyear 10000)) 100)))
	 (jdday (- sjd (+ (* jdyear 10000) (* jdmo 100))))
	 (modays (nth (1- jdmo) *jdmodays*))
	 )
    (incf jdday)
    (when (> jdday modays)
      (setq jdday 1)
      (incf jdmo))
    (when (> jdmo 12)
      (setq jdmo 1)
      (incf jdyear)
      )
    (read-from-string
     (format nil "~4d~2,'0d~2,'0d" jdyear jdmo jdday))
    ))

(defun difference-jeff-dates (jd1 jd2)
  (- (norm-jeff-date jd1) (norm-jeff-date jd2)))

;;; Like APL compress, uses a list of t/nil to select from a target of
;;; the same length.  This is a really inefficient implementation of
;;; this.  Ought to be destructive some day.

(defun compress (selector target &aux r)
  (mapcar #'(lambda (a b) (if a (push b r))) selector target)
  (reverse r))

;;; Why isn't this standard in lisp???

(defun flatten (l)
  (cond ((null l) ())
	((atom l) (list l))
	(t (append (flatten (car l))
		   (flatten (cdr l))))))

;;; Copy the first c lines of a file (skipping s lines first).  Used
;;; to break up huge files.

(defun fscopy (infile outfile c &optional (s 0))
  (with-open-file (instream infile :direction :input)
    (with-open-file (outstream outfile :direction :output :if-exists :supersede)
      (when (not (zerop s))
	(loop for i from 1 to (* s 2) ; count the line-turn
	      do (read-line instream))
	)
      (loop for i from 1 to c
	    do 
	    (princ (read-line instream) outstream)
	    (read-line instream)
	    (terpri outstream)
	    )
      )))

;;; ===================================================================
;;; Outdated or un-understood stuff.

#| Some versions of lisp require ^M processing, so I've dyked this, 
   and I can read it into individual programs and modified as needed.  

(defun scan-stream-for-line-beginning (stream string &key (skip 0))
  (loop with string-length = (length string)
        with sl+ = (+ skip string-length)
	for line = (read-line! stream)
	until (or (null line)
		  (and (>= (length line) sl+)
		       (string-equal (subseq line skip sl+) string)
		       ))
	finally (return line)))

(defun scan-stream-for-line-beginning (stream string &key (skip 1))
  (loop with string-length = (length string)
        with sl+ = (+ skip string-length)
	for line = (read-line stream nil nil)
	until (or (null line)
		  (and (>= (length line) sl+)
		       (string-equal (subseq line skip sl+) string)
		       ))
	finally (return line)))

|#

(defun read-line! (s)
  (loop with c* 
	for c = (read-char s nil nil)
	do
	(cond ((null c) (return nil))
	      ((member c '(#\linefeed #\newline #\return) :test #'char-equal)
	       (return (coerce (reverse c*) 'string)))
	      (t (push c c*)))))

;;; This writes a hashtable to a file

(defun write-hashtable (file ht)
  (with-open-file (stream file :direction :output :if-exists :supersede)
		  (maphash #'(lambda (key value)
			       (print (cons key value) stream)) 
			   ht)))

;;; This reads a hashtable from a file

(defun read-hashtable (file)
  (let* ((ht (make-hash-table :test #'equal)))
    (with-open-file (stream file :direction :input)
		    (loop for line = (read stream nil nil)
			  until (null line)
			  do
			  (push (cdr line) (gethash (car line) ht))))
    ht))

;;; ===================================================================
#| 

20060924: think that this is massively out-of-date, and superseded by
apply-across-dirtree.

;;; Gets the names from the current directory through ls and parses
;;; them into a tree which can be displayed in menu dialog format in
;;; order to choose from among the files.  The only function that this
;;; file provides is histogramming the results.

(defun test ()
  (load-and-plot (select-file "/data/96mar27/collected")))

;;; Probably you'll want to do the parts of this separately in your
;;; own code, but it's here as example.  Note that the dir is for
;;; filer so you have to pass it without the final '/'

(defun select-file (dir)
  (get-dir-filenames dir) ; saves out into *d*
  (treeify (mapcar #'parse-filename *d*)) ; saves out into *dirtree*
  (selfn *dirtree*)
  )

;;; This uses filer to put quotes around the filenames.  It's a wee
;;; bit obscure, and won't work for anyone who doesn't have filer.
;;; The better way to do it would be to create the dir listing and
;;; then read it in line by line, but that takes too long.  The
;;; filenames end up looking like this: " path/namepart.namepart " and
;;; have to be trimmed.

(defvar *d* ()) ; gets the list of filenames as strings

(defun get-dir-filenames (dir)
  (system (format nil "filer -d ~a -c \"\\\"\" -m \"*\" -r \"\\\"\" > /tmp/foo" dir))
  (setq *d* (car (read-data-columns "/tmp/foo" 1)))
  (length *d*))

;;; Parse up filesnames at dots.We get back a list like: (parsedfn...)
;;; where parsedfn is: ("fullname" "namepart1" "namepart2"...)
;;; The tree gets built from these later.

(defun parse-filename (fn)
  (pf2 (string-trim ". " fn))) ; removes filer leftovers . and spaces

(defun pf2 (fn)
  ;; the returned form is (fullnamestring (firstpart secondpart ...))
  (list fn (reverse (pf3 fn (length fn) "" ()))))

(defun pf3 (fn lfn name store)
  (cond ((zerop lfn) (push name store)); end of the line
	((char= #\. (aref fn 0)); looking at a dot, nextname, please
	 (pf3 (string-left-trim "." fn) (decf lfn) "" (push name store)))
	(t ; okay, not done and note a ., recurse with this name
	 (pf3 (subseq fn 1)
	      (decf lfn)
	      (format nil "~a~a" name (aref fn 0))
	      store))
	))
	 
;;; Take the list of names and their parses produced by parse-filename
;;; (applied to the list of names) and produce the full tree from it.
;;; The tree is stored in *dirtree* as an alist.

(defvar *dirtree* ())

(defun treeify (fnl)
  (setq *dirtree* ())
  (mapcar #'tfy2a fnl))

;;; Insert each filename (as parsed) into the tree.

(defun tfy2a (fn)
  (tfy2b (cadr fn) (car fn)))

(defun tfy2b (fnparts fullname)
  (setq *dirtree* (tfy3 (car fnparts) (cdr fnparts) *dirtree* fullname)))

;;; Here's the hard part.  We are handed the current part of the filename,
;;; the rest of the filename, the level of the tree that this current part
;;; should go into, and the full filename.

(defun tfy3 (fnpart fnrest dtlevel fullname)
  (cond ((null fnpart) ; If we're through...
	 (cons (list "Done" fullname) dtlevel)) ; ... return the full filename
				  ; jammed into the front of 
                                  ; the current level list.
        ; Otherwise, find this part in the current level alist.
	(t (let ((loc (assoc fnpart dtlevel :test #'string=)))
	     (cond ((null loc) ; it's not there yet in this level
		    ;; make a new entry and recurse
		    (cons (cons fnpart
                               (tfy3 (car fnrest) (cdr fnrest) () fullname))
			  dtlevel))
		   (t ; else, it is there, just recursewith the current loc
;		    (print loc)
		    (cons (cons (car loc)
				(tfy3 (car fnrest) (cdr fnrest) (cdr loc) fullname))
			  (remove loc dtlevel))
		   ))))
	))


;;; Select a filename from the dirtree. 

(defun selfn (tr)
  (cond ((stringp (car tr)) (car tr))
	(t (selfn (cdr (nth (seln (mapcar #'car tr)) tr))))))

(defun seln (atoms)
  (choose-item-dialog "Select a Condition Part" atoms))

;;; Apply himage to the chosen file and plot it.

(defvar *histdata* ())

(defun load-and-plot (fn)
  (system (format nil "/users/hahn/bin/himage ~a > /tmp/histdata" fn))
  (setq *histdata* (read-data-columns "/tmp/histdata" 2))
  ;; remove 0 (first) entry because it's the 0phase, which has all the black points!!
  (pop (car *histdata*))
  (rplacd *histdata* (list (cdadr *histdata*)))
  (plot-points (car *histdata*)(cadr *histdata*)) 
  )

;;; Another version of the same thing, also massively out of date!

;;; --- The fn (understand-directory "path") returns a list of file
;;; structures for your computing pleasure.  A file structure gets
;;; created for each entry.  Most of the stuff in the long ls listing
;;; gets ignored in this process.

(defstruct file size month day time name directory?)

;;; Get the contents of a directory.  This returns a list which is the result
;;; of ls -l on the given path.  

;;; (God Damned F--king Unix is inconsistent about whether the total
;;;  size is given or not in ls -l.  When you give it a path that's a
;;;  directory it gives you the total, but if you don't, it doesn't.
;;;  So I have to see if the first thing is a total line, and either
;;;  remove it or not as the case.  Foo!)

(defun understand-directory (path &aux result)
  ;; could mapcar this but we have to check for each case whether it 
  ;; conforms to what we're expecting, that is, contains 9 items, else
  ;; ignore it.  (The "total" line in ls has only 2, e.g.,)
  (dolist (item (ud2 path))
    (if (= 9 (length item))
	(push (make-file-struct item) result))
    ;; look for links
    (if (and (= 11 (length item)) (string-equal "->" (nth 9 item)))
	(push (make-file-struct item) result))
    )
  result)

(defun ud2 (path &aux allfiles tmpfile)
    (setq tmpfile (format nil "/tmp/~a.dirinfo.~a" 
			  (system::getenv "LOGNAME")
			  (random 10000)))
    (system (format nil "ls -l ~a > ~a" path tmpfile))
    (with-open-file (f tmpfile :direction :input)
      (prog (onefile item)
	loop
	    (setq onefile nil)
	    (setq item (read-line f nil nil))
	    (if (null item) (return t))
	    (dolist (item (parse-string item :make-substitutions nil))
	      (push item onefile))
	    (push (reverse onefile) allfiles)
	    (go loop)
	    ))
    (delete-file tmpfile)
    allfiles
    )

(defun make-file-struct (flist)
  (make-file
             :size (read-from-string (nth 4 flist))
	     :name (nth 8 flist)
	     :month (read-from-string (nth 5 flist))
	     :day (read-from-string (nth 6 flist))
	     ;; Time could be the year in older files.  In this case,
	     ;; parse-time will return a LIST instead of a numberp so that
	     ;; other can tell the difference.
	     :time (parse-time (nth 7 flist))
	     :directory? (if (equal #\d (aref (nth 0 flist) 0)) t nil)
	     ))

;;; compute the number of minutes after midnight.  Note that if it's
;;; an older file, the time will be a year, in which case we simply
;;; return it as a LIST of the time, which gives

(defun parse-time (time)
  (if (= 4 (length time)) ; it's like 1995...
      (list (read-from-string time)) ;... so just return that as a list.
   (let* ((string (format nil "~a" time))
	  (hours (read-from-string (subseq string 0 2)))
	  (mins (read-from-string (subseq string 3 5)))
	  )
     (+ (* hours 60) mins))))

;;; Some utilities.

(defun sort-files-by-time (files)
  (sort files #'(lambda (a b) (< (file-time a) (file-time b)))))

;;; prob. not needed unless not in commonlisp -- remove well after 961001
;;; This is used by any function that is going to do read-line's and
;;; needs to be called both interactively and sometimes from another
;;; function.  The problem is that when the user type the call at the
;;; lisp reader, there's a newline in the input buffer, but when it's
;;; called from some other function that's already doing reads, there
;;; isn't.  This missynchs input.  This function eats a single newline
;;; that's left in the input buffer, otherwise does nothing.

(defun clear-input ()
  (if (equal (peek-char) #\newline)
      (read-line)))

|#

;;; ===================================================================
;;; Takes a connectivity matrix (here represented as a list of:
;;; (... (from-node to-node-1 to-node-2 ...) ...) where the nodes are
;;; just integers (although this shouldn't matter, and produces a new
;;; list of: (...(new-node-number member-node-1 member-node-2...) ...)
;;; where all the nodes that connect to one another are grouped.

(defvar *test-matrix* nil)
(defvar *colors->contigs* (make-hash-table))
(defvar *contigs->colors* (make-hash-table))

(defun unify-dups (connections &aux (current-color 0))
  (clrhash *colors->contigs*)
  (clrhash *contigs->colors*)
  (loop with report-every = (round (/ (length connections) 10))
	for (from-contig . to-contigs) in connections
	as from-already-colored = (gethash from-contig *contigs->colors*)
	as k from 1 by 1
	do (when (zerop (mod k report-every)) (format t "@~a color=~a~%" k current-color))
	(cond (from-already-colored
		  (loop for to-contig in to-contigs
			do (recolor-all to-contig from-already-colored)))
		 (t (let ((from-already-colored (incf current-color)))
		      (setf (gethash from-contig *contigs->colors*) from-already-colored)
		      (setf (gethash from-already-colored *colors->contigs*) (list from-contig))
		      (loop for to-contig in to-contigs ; at the moment we assume that all are "perfect" links
			    do (recolor-all to-contig from-already-colored)
			    (pushnew to-contig (gethash from-already-colored *colors->contigs*) :test #'equal)
			    )
		      )))
	)
  (rebuild-colors->contigs-table)
  )

;;; The colors->contigs table is highly redundant, and is rebuilt from 
;;; the opposite table so that it can be properly understod AFTER recoloring.

(defun rebuild-colors->contigs-table ()
  (clrhash *colors->contigs*)
  (maphash #'(lambda (contig color)
	       (pushnew contig (gethash color *colors->contigs*) :test #'equal))
	   *contigs->colors*))

(defun recolor-all (contig color)
  (let ((current-color (gethash contig *contigs->colors*)))
    (cond ((and current-color (= color current-color)) nil)
	  (t (setf (gethash contig *contigs->colors*) color)
	     (pushnew contig (gethash color *colors->contigs*))
	     (loop for link-contig in (gethash current-color *colors->contigs*)
		   do (recolor-all link-contig color))))))

(defun test-unify-dups (&key (outfile "c:/temp/uni.out") (n-nodes 1000) (xns-per-node 2))
  (gen-test-xm n-nodes xns-per-node)
  (unify-dups *test-matrix*)
  ;; Reporting
  (with-open-file (o outfile :direction :output :if-exists :supersede)
    (loop for color being the hash-keys of  *colors->contigs*
	  using (hash-value contigs)
	  as k from 1 by 1
	  do (format o "~a = ~a~%" color contigs)
	  finally (format t "There where ~a groups.~%" k))))

(defun gen-test-xm (n-nodes xns-per-node)
  (setq *test-matrix*
	(loop for node from 1 to n-nodes
	      collect (loop for xn from 1 to xns-per-node
			    collect (random n-nodes)))))

(defun categorize (objects function &key (count? nil) (FAST-BUT-REDUNDANT? t))
  "Given a list of object, and a function that applies to those objects and returns
 a single category or list of categories re-organize the objects by the cateogies.
 This uses EQUAL hash tables, so the function must return something that can be EQUAL tested.
 Example: (categorize '(1 2 3 4 5) #'evenp) => ((NIL 5 3 1) (T 4 2))
 [The order of the result categorize and objects in each is indeterminate.]
 More complex example:
 (categorize '(\"Jeff\" \"Joe\" \"Sam\" \"Sall\") #'(lambda (object) (list (elt object 0) (length object))))
 ((3 \"Sam\" \"Joe\") (4 \"Sall\" \"Jeff\") (#\J \"Joe\" \"Jeff\") (#\S \"Sall\" \"Sam\"))
 Note that we get BOTH categories in this case by first letter AND by length.
 By default this uses PUSH to add object to the category, which is fast,
 but doesn't protect against redundant entries. If you want to protect for this,
 add the keyword :FAST-BUT-REDUNDANT? NIL (default is T)
 By default you get all the entries back, but sometimes you just want a count.
 Adding the keyword :COUNT? t (default = nil) will spare you the pain of the whole list!"
  (let ((table (make-hash-table :test #'equal)))
    (loop for object in objects
	  as cats = (let ((cats (funcall function object)))
		      (if (and (listp cats) (not (null cats)))
			  cats (list cats)))
	  do (loop for cat in cats
		   do (if FAST-BUT-REDUNDANT? 
			  (push object (gethash cat table))
			(pushnew object (gethash cat table) :test #'equal))))
    (loop for k being the hash-keys of table
	  using (hash-value v)
	  collect 
	  (if count?
	      (list k (length v))
	    (cons k v)))))
    
(defun string-sub (from to in)
  "Substitute to for from every time it appears in in. CAREFUL THAT YOU DON'T CREATE AN INFINITE SUBSTITUTION LOOP!"
  (let ((p (search from in)))
    (if p
        (string-sub from to
                    (format nil "~a~a~a" (subseq in 0 p) to (subseq in (+ p (length from)))))
      in)))