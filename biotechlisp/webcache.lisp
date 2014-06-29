;;; Can be used in ACL in place of CG:WEB-PAGE-CONTENTS to look first in our
;;; local cache (indicated by cache-dir) and then whereever.  When the page
;;; is loaded, it gets cached.

(defvar *cache-dir->cache-table-table* (make-hash-table :test #'equal))

(defun WEB-PAGE-CONTENTS-CACHED (url cache-dir &key (reload? nil) (trace? t))
  (cond ((or reload? (null (get-cached-page url cache-dir)))
	 (when trace? (format t "Going to the internet for ~s~%" url))
	 (cache-web-page url (CG:WEB-PAGE-CONTENTS url) cache-dir))
	(t (when trace? (format t "Going to local cache for ~s~%" url))
	   (get-cached-page url cache-dir))))

(defun get-cached-page (url dir)
  (let ((cache-table (ensure-cache-table dir)))
    (load-from-cache-file (gethash url cache-table) dir)))

(defun load-from-cache-file (filename dir)
  (let ((full-filename (format nil "~a/~a.txt" dir filename)))
    (when (probe-file full-filename)
	  (with-open-file (i full-filename)
			  (read i)))))

(defun ensure-cache-table (dir)
  (or (gethash dir *cache-dir->cache-table-table*)
      (setf (gethash dir *cache-dir->cache-table-table*)
	    (load-or-create-cache-index dir))))

(defun load-or-create-cache-index (dir)
  (let ((index-filename (form-index-filename dir)))
    (cond ((probe-file index-filename) (load-index index-filename))
	  (t (create-index index-filename)
	     (load-index index-filename)))))

(defun form-index-filename (dir)
  (format nil "~a/cache.index" dir))

;;; And index entry is ("url" . idnumber)

(defun load-index (index-filename)
  (with-open-file (i index-filename)
    (loop with index-table = (make-hash-table :test #'equal)
	  as entry = (read i nil nil)
	  until (null entry)
	  do (setf (gethash (car entry) index-table) (cdr entry))
	  finally (return index-table))))

(defun create-index (index-filename)
  (with-open-file (o index-filename :direction :output)
		  ))
    
;;; Add a new entry into the cache.  The id is just the hash table count,
;;; and the page just gets stored in the cache-dir under that number, and a
;;; line added to the index.  If the page is already in the table, then 
;;; this is a reload, and we reuse its number.

(defun cache-web-page (url contents cache-dir)
  (let* ((table (ensure-cache-table cache-dir))
	 (existing-id (gethash url table)) 
	 (next-id (1+ (hash-table-count table))))
    (cond (existing-id (store-cache-page existing-id contents cache-dir))
	  (t (store-cache-page next-id contents cache-dir)
	     (with-open-file (o (form-index-filename cache-dir) :direction :output :if-exists :append)
			     (print (cons url next-id) o))
	     (setf (gethash url table) next-id))
	  ))
  contents
  )

(defun store-cache-page (id contents dir)
  (with-open-file (o (format nil "~a/~a.txt" dir id) :direction :output)
    (print contents o)))
  