; This file contains additional Lisp utilities that were not included 
; in _On Lisp_ or _ANSI Common Lisp_.

; This code is copyright 1995 by Paul Graham, but anyone who wants to
; use it is free to do so.  Report bugs to lispcode@paulgraham.com.


; Lists

(defun agree (x y)
  (if (or (null x) (null y))
      t
      (if (equal (car x) (car y))
          (agree (cdr x) (cdr y)))))

(defun assocify (source)
  (labels ((rec (source acc)
             (let ((rest (cddr source)))
               (if (consp rest)
                   (rec rest (cons (cons (car source) (cadr source)) acc))
                   (nreverse (cons (cons (car source) (cadr source)) acc))))))
    (if source (rec source nil) nil)))

(defun duplicate (obj lst &key (test #'eql))
  (member obj (cdr (member obj lst :test test))
          :test test))

(defun careql (x y)
  (and (consp x) (eql (car x) y)))

(defmacro carin (x &rest args)
  (with-gensyms (g)
    `(let ((,g ,x))
       (and (consp ,g) (in (car ,g) ,@args)))))

(defun carat (x)
  (if (consp x) (car x) x))

(define-modify-macro pushend (obj)
  (lambda (place obj)
    (nconc place (list obj))))

(define-modify-macro merge-into (obj fn)
  (lambda (place obj fn)
    (merge 'list place (list obj) fn)))

(defun firstn (lst n)
  (if (or (null lst) (<= n 0))
      nil
      (cons (car lst)
            (firstn (cdr lst) (- n 1)))))

(defun butlastn (seq n)
  (subseq seq 0 (- (length seq) n)))

(defun rotlist (x)      
  (if (cdr x) 
      (cons (car (last x)) (butlast x))
      x))

(defun prefix (pref str)
  (search pref str :end2 (min (length pref) (length str))))
    
(defun suffix (str suff)
  (search suff str :start2 (- (length str) (length suff))))

(defun insert-before (before after lst)
  (cond ((null lst) nil)
        ((eql (car lst) after)
         (cons before lst))
        (t (cons (car lst) (insert-before before after (cdr lst))))))

(defun insert-after (before after lst)
  (cond ((null lst) nil)
        ((eql (car lst) before)
         (cons before (cons after (cdr lst))))
        (t (cons (car lst) (insert-after before after (cdr lst))))))

(defmacro pull-nth (n place)
  (multiple-value-bind (vars forms var set access)
                       (get-setf-method place)
    (let ((g (gensym)))
      `(let* ((,g ,n)  
              ,@(mapcar #'list vars forms)
              (,(car var) (delete-nth ,g ,access)))
         ,set))))
     
(defun delete-nth (n lst)
  (cond ((< n 0) (error "Bad arg to delete-nth"))
        ((= n 0) (cdr lst))
        (t (let ((rest (nthcdr (1- n) lst)))
             (pop (cdr rest))
             lst))))

(defmacro push-nth (n obj place)
  (multiple-value-bind (vars forms var set access)
                       (get-setf-method place)
    (with-gensyms (g h) 
      `(let* ((,g ,n)
              (,h ,obj)
              ,@(mapcar #'list vars forms)
              (,(car var) (ninsert-nth ,g ,h ,access)))
         ,set))))
    
(defun ninsert-nth (n obj lst)
  (if (< n 0)
      (error "Bad arg to ninsert-nth")
      (let ((rest (nthcdr n lst)))
        (push obj (cdr rest))
        lst)))
  
(defun insert-elt-after (elt ins lst)
  (if (null lst)
      nil     
      (if (eql (car lst) elt)
          (cons (car lst) (cons ins (cdr lst)))
          (cons (car lst) (insert-elt-after elt ins (cdr lst))))))


; Control

(defmacro bind (&rest args)
  `(multiple-value-bind ,@args))

(defun ifnot (bad val) 
   (unless (eql bad val) val))

(defmacro nullor (x y)
  (with-gensyms (g)
    `(let ((,g ,x))
       (if (zerop (length ,g)) ,y ,g))))

(defmacro until (test &rest body)
  `(do ()
       (,test)
     ,@body))

(defmacro casequal (val &rest clauses)
  (let ((g (gensym)))
    `(let ((,g ,val))
       (cond ,@(mapcar #'(lambda (cl)
                           `(,(if (eql (car cl) t)
                                  t
                                  `(equal ,g ,(car cl)))
                             ,@(cdr cl)))
                       clauses)))))


; Iteration

(defmacro do-all (var x &rest body)
  (with-gensyms (g)
    `(let ((,g ,x))
       (if (consp ,g)
           (dolist (,var ,g) ,@body)
           (let ((,var ,g)) ,@body)))))

(defmacro dolists (pairs &rest body)
  (with-gensyms (f)
    (let ((parms (mapcar #'(lambda (p) (gensym)) pairs)))
      `(labels ((,f ,parms
                  (when (or ,@parms)
                    (let ,(mapcar #'(lambda (p g)
                                      (list (car p) `(car ,g)))
                                  pairs
                                  parms)
                      ,@body)
                    (,f ,@(mapcar #'(lambda (g) `(cdr ,g))
                                  parms)))))
         (,f ,@(mapcar #'cadr pairs))))))

(defmacro do3 (v1 v2 v3 list &rest body)
  (with-gensyms (g h)
    `(let ((,g ,list))
       (do ((,h ,g (cdr ,h)))
           ((null ,h) nil)
         (let ((,v1 (car ,h))
               (,v2 (if (cdr ,h) (cadr ,h) (car ,g)))
               (,v3 (if (cdr ,h)
                        (if (cddr ,h)
                            (third ,h)
                            (car ,g))
                        (if (cdr ,g)
                            (second ,g)
                            (car ,g)))))
           ,@body)))))

; Assumes 3 args.  Inefficient.

(defmacro do-cyclic (parms source &rest body)
  (let ((s (gensym)))
    `(let ((,s ,source))
       (case (length ,s)
         (0 nil)
         (1 (let ((,(first parms) (first ,s))
                  (,(second parms) (first ,s))
                  (,(third parms) (first ,s)))
              ,@body))
         (2 (let ((,(first parms) (second ,s))
                  (,(second parms) (first ,s))
                  (,(third parms) (second ,s)))
              ,@body)
            (let ((,(first parms) (first ,s))
                  (,(second parms) (second ,s))
                  (,(third parms) (first ,s)))
              ,@body))
         (t (do-tuples/c ,parms (rotlist ,s)
              ,@body))))))

(defmacro do-plist (v1 v2 plist &rest body)
  (with-gensyms (rec rest pl)
    `(labels ((,rec (,v1 ,v2 ,rest)
                 ,@body
                 (when ,rest
                   (,rec (car ,rest) (cadr ,rest) (cddr ,rest)))))
       (let ((,pl ,plist))
         (when (consp ,pl)
           (,rec (car ,pl) (cadr ,pl) (cddr ,pl)))))))


; Files & Streams

(defconstant eof (gensym))

(defun clear-file (name)
  (aif (probe-file name) (delete-file it)))

(defmacro do-chars (var str &rest body)
  (with-gensyms (g) 
    `(let ((,g ,str))
       (do ((,var (read-char ,g nil eof) (read-char ,g nil eof)))
           ((eql ,var eof)) 
         ,@body))))

(defmacro do-stream (var str &rest body)
  (with-gensyms (g)
    `(let ((,g ,str))
       (do ((,var (read ,g nil eof) (read ,g nil eof)))
           ((eql ,var eof))
         ,@body))))

(defun file-read (path)
  (awhen (probe-file path)
    (with-infile str it
      (read str nil nil))))
             
(defun file-read-line (path)
  (awhen (probe-file path)  
    (with-infile str it     
      (read-line str nil nil))))
                        
(defun file-contents (path) 
  (awhen (probe-file path)
    (with-infile str it
      (let ((sb (make-sb))) 
        (do-chars c str
          (sb-append sb c))
        sb))))

(defun file-to-stream (path &optional (out *standard-output*))
  (awhen (probe-file path)      
    (with-infile str it
      (do-chars c str
        (princ c out)))))
       
(defun file-rest (path)
  (awhen (probe-file path)
    (with-open-file (str it :direction :input)
      (read-rest str))))

(defun read-rest (str)
  (let ((x (read str nil eof)))
    (if (eql x eof)
        nil
        (cons x (read-rest str)))))

(defun suck-dry (str)
  (do ((c (read-char str nil :eof) (read-char str nil :eof)))
      ((eql c :eof) nil)))

(defmacro do-lines (var path &rest body)
  (with-gensyms (p str)
    `(let ((,p (probe-file ,path)))
       (when ,p
         (with-infile ,str ,p
           (do ((,var (read-line ,str nil eof)
                      (read-line ,str nil eof)))
               ((eql ,var eof))
             ,@body))))))

(defun file-lines (path)
  (let ((count 0))
    (do-lines line path
      (incf count))
    count))

(defmacro with-infile (var fname &rest body)
  (with-gensyms (v c f)
    `(let ((,f ,fname))
       (in-case-error
         (with-open-file (,var ,f :direction :input)
           ,@body)
         (format *error-output* "Error reading from ~s.~%" ,f)))))

(defmacro with-outfile (var fname &rest body)
  (with-gensyms (v c f)
    `(let ((,f ,fname))
       (in-case-error
         (with-open-file (,var ,f :direction :output
                                  :if-exists :supersede)
           ,@body)
         (format *error-output* "Error writing to ~s.~%" ,f)))))

(defmacro with-outfiles (pairs &rest body)
  (if (null pairs)
      `(progn ,@body)
      (if (oddp (length pairs))
          (error "Odd length arg to with-outfiles")
          `(with-outfile ,(car pairs) ,(cadr pairs)
             (with-outfiles ,(cddr pairs) ,@body)))))

(defun copy-file (from to)
  (with-open-file (in from :direction :input
                           :element-type 'unsigned-byte)
    (in-case-error
      (with-open-file (out to :direction :output
                              :element-type 'unsigned-byte)
        (do ((i (read-byte in nil -1)
                (read-byte in nil -1)))
            ((minusp i))
          (declare (fixnum i))
          (write-byte i out)))
      (format *error-output* "Error writing to ~s.~%" to))))


; Strings

(defun make-sb ()
  (make-array '(0) :element-type 'string-char :adjustable t
          :fill-pointer 0))     
    
(defmacro sb-append (as c)
  `(vector-push-extend (the character ,c) (the string ,as)))

(defun read-as-string (in)
  (with-output-to-string (out)
    (do ((c (read-char in nil :eof) (read-char in nil :eof)))
        ((eql c :eof))
      (write-char c out))))

(defun read-list-from-string (str &optional (start 0))
  (bind (val pos) (read-from-string str nil eof :start start)
    (if (eql val eof)
        nil
        (cons val (read-list-from-string str pos)))))

(defun read-numbers-from-string (str &optional (start 0))
  (bind (val pos) (ignore-errors
                    (read-from-string str nil eof :start start))
    (if (numberp val)
        (cons val (read-numbers-from-string str pos))
        nil)))

(defun read-number-from-string (string)
  (let ((n (ignore-errors (read-from-string string))))
    (if (numberp n) n nil)))

(defun string->hex (str)
  (with-output-to-string (out)
    (dotimes (i (length str))
      (let ((c (char str i)))
        (format out "~2,'0X" (char-code c))))))
              
(defun hex->string (str)
  (if (evenp (length str))
      (with-output-to-string (out) 
        (let ((*read-base* 16))
          (dotimes (i (/ (length str) 2)) 
            (princ (code-char (read-from-string str  nil nil
                                                :start (* i 2)
                                                :end (+ (* i 2) 2)))
                   out))))
      (error (format nil "odd-length hex string: ~A" str))))


; Hash Tables

(defun nonempty-ht (ht)
  (maphash #'(lambda (k v) (return-from nonempty-ht t))
           ht)
  nil)

(defun ht-car (ht)
  (maphash #'(lambda (k v) (return-from ht-car v))
           ht))

(defun hash-keys (ht)
  (let ((acc nil))
    (maphash #'(lambda (k v) 
                 (declare (ignore v)) 
                 (push k acc))
             ht)
    acc))

(defun hash-vals (ht)
  (let ((acc nil))
    (maphash #'(lambda (k v)
                 (declare (ignore k))
                 (push v acc))
             ht)
    acc))

(defun hash-pairs (ht)
  (let ((acc nil))
    (maphash #'(lambda (k v)
                 (push (cons k v) acc))
             ht)
    acc))

(defun somehash (fn ht)
  (maphash #'(lambda (k v)
               (when (funcall fn v)
                 (return-from somehash v)))
           ht)
  nil)
     
(defun key-match (ht1 ht2)
  (maphash #'(lambda (k v)
               (declare (ignore v))
               (when (gethash k ht2)
                 (return-from key-match k)))
           ht1)
  nil)

(defun write-hash (ht str)
  (maphash #'(lambda (k v)
               (print (cons k v) str))
           ht))

(defun read-hash (ht str)
  (do-stream pair str
    (setf (gethash (car pair) ht)
          (cdr pair)))
  ht)


; Errors & Debugging

(defun ero (&rest args)
  (print (if (cdr args) args (car args))
         *error-output*))

(defmacro safely (expr)
  (with-gensyms (ret err)
    `(bind (,ret ,err) (ignore-errors ,expr)
       (if (typep ,err 'error)
           nil
           (values ,ret ,err)))))

(defmacro in-case-error (expr err)
  (with-gensyms (val cond)
    `(bind (,val ,cond) (ignore-errors ,expr)
       (if (typep ,cond 'error)
           (progn
             ,err
             (error ,cond))
           ,val))))


; Date & Time

(defun time-string ()
  (multiple-value-bind (s m h) (get-decoded-time)
    (format nil "~A:~2,,,'0@A:~2,,,'0@A" h m s)))

(defconstant months
  '#("Jan" "Feb" "Mar" "Apr" "May" "June" "Jul" "Aug" "Sept" "Oct" "Nov" "Dec"))

(defun date-string ()
  (multiple-value-bind (ig no re d mo y) (get-decoded-time)
    (declare (ignore ig no re))
    (format nil "~A ~A ~A" d (svref months (1- mo)) y)))

(defun date+time-string (&optional (u (get-universal-time)))
  (multiple-value-bind (s m h d mo y) (decode-universal-time u)
    (format nil "~A ~A ~A ~A:~2,,,'0@A:~2,,,'0@A"
                d (svref months (1- mo)) y h m s)))


; String Processing

(defun white? (c)
  (in c #\Space #\Tab #\Return #\Newline))

(defun constituent (c)
  (and (graphic-char-p c)
       (not (char= c #\  ))))

(defconstant whitechars '#(#\Space #\Tab #\Return #\Newline))

(defun extract-paragraphs (string &optional (pos 0))
  (if (= pos (length string))
      nil
      (let ((n (blank-line string pos)))
        (if n
            (cons (nonwhite-substring string pos n)
                  (extract-paragraphs string (1+ n)))
            (if (find-if-not #'white? string :start pos)
                (list (nonwhite-substring string pos (length string)))
                nil)))))

(defun nonwhite-substring (string start end)
  (if (find-if-not #'white? string :start start :end end)
      (progn
        (while (white? (char string start))
          (incf start))
        (while (white? (char string (- end 1)))
          (decf end))
        (subseq string start end))
      ""))

(defun blank-line (string &optional (start 0))
  (let ((n (linebreak-pos string start)))
    (if (and n (not (= n (1- (length string)))))
        (let ((n2 (linebreak-pos string (1+ n))))
          (if n2
              (if (find-if-not #'white? string :start n :end n2)
                  (blank-line string n2)
                  n2)
              nil))
        nil)))

; A linebreak is either a linefeed, or a cr-linefeed.  Returns
; the position of the last char in the next linebreak after start.

(defun linebreak-pos (string &optional (start 0))
  (or (position #\Newline string :start start)
      (let ((n (position #\Return string :start start)))
        (when n
          (if (and (not (= n (1- (length string))))
                   (eql (char string (1+ n)) #\Newline))
              (1+ n)
              n)))))

(defun extract-lines (string &optional (pos 0))
  (if (= pos (length string))
      nil
      (let ((n (linebreak-pos string pos)))
        (if n
            (cons (nonwhite-substring string pos n)
                  (extract-lines string (1+ n)))
            (list (nonwhite-substring string pos (length string)))))))


(defun extract-tokens-if (str sep-test &optional (start 0))
  (let ((p1 (position-if-not sep-test str :start start)))
   (if p1
       (let ((p2 (position-if sep-test str :start p1)))
         (cons (subseq str p1 p2)
               (if p2
                   (extract-tokens-if str sep-test p2)
                   nil)))
       nil)))

(defun separated-tokens (str sep)
  (mapcar #'(lambda (tok)
              (string-trim whitechars tok))
          (extract-tokens-if str
                             #'(lambda (c) (eql c sep)))))


(defun extract-tokens (str &optional (start 0))
  (let ((p1 (position-if #'constituent str :start start)))
   (if p1
       (if (eql (char str p1) #\")
           (if (< p1 (- (length str) 1))
               (let ((p2 (position #\" str :start (1+ p1))))
                 (if (and p2 (< p2 (- (length str) 1)))
                     (cons (string-trim "\"" (subseq str p1 (1+ p2)))
                           (extract-tokens str (1+ p2)))
                     (list (string-trim "\"" (subseq str p1)))))
               nil)
           (let ((p2 (position-if-not #'constituent
                                      str :start p1)))
             (cons (subseq str p1 p2)
                   (if p2
                       (extract-tokens str p2)
                       nil))))
       nil)))

(defun first-token (str)
  (let ((p1 (position-if #'constituent str)))
    (if p1
        (subseq str p1 (position-if-not #'constituent
                                        str :start p1))
        nil)))
