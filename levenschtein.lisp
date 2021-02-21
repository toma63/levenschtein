
(load "../heap-cl/heap.lisp")

(defclass partial ()
  ((moves :initform '()
	  :accessor moves
	  :initarg :moves)
   (cost :initform 0 :accessor cost :initarg :cost)
   (txt-pos :initform 0 :accessor txt-pos :initarg :txt-pos)
   (max-txt :initform 0 :accessor max-txt :initarg :max-txt)
   (pat-pos :initform  0 :accessor pat-pos :initarg :pat-pos)
   (max-pat :initform 0 :accessor max-pat :initarg :max-pat)
   (matched :initform 0 :accessor matched :initarg :matched) ; number of matched positions
   (distance :initform 0 :accessor distance :initarg :distance)
   (delta :initform 0 :accessor delta :initarg :delta))
  (:documentation "a partial solution to a string edit distance"))

(defmethod print-object ((obj partial) stream)
      (print-unreadable-object (obj stream :type t)
        (with-accessors ((moves moves)
                         (cost cost)
			 (txt-pos txt-pos)
			 (pat-pos pat-pos)
			 (matched matched)
			 (distance distance)
			 (delta delta))
            obj
          (format stream "moves: ~a~% cost: ~a~% matched: ~a~% txt-pos: ~a~% pat-pos: ~a~% distance: ~a~% delta: ~a~%"
		  moves cost matched txt-pos pat-pos distance delta))))

;; simple test data
(defvar txt "abcdefghij")
(defvar pat "abdqfghijk")
(defvar thou "thou shalt not")
(defvar you "you should not")

;; related DNA sequences
(defvar tseq "ATCGTATCAGAGCTTTCGTCTGATGCTCGATTGTAA")
(defvar tal1    "GTATCAGAGCTTTCGTCTGATGCTCGATTGTAA")
(defvar tal2 "ATCGTATCTTTCGTCTGATGCTCGATTGTAA")
(defvar tal3 "ATCGTAGAGCTTTCGTCTGATGCTCGATTGTAA")
(defvar tal4 "ATCGTATCAGAGCGTCTGATGCATTGTAA")
(defvar tal5 "ATCGTAGAGCCGTCTGCTCGATTGTAA")
(defvar tal6 "GCTTTCGTCTGATGCTCGATTGTAACAGCAGCAGCAG")

(defgeneric compare (obj1 obj2)
  (:documentation "comparison function for use by a heap"))

;; sorting compare function uses summes cost, distance, delta less match distance
(defmethod compare ((obj1 partial) (obj2 partial))
  (with-slots ((cost1 cost) (matched1 matched) (delta1 delta) (distance1 distance)) obj1
    (with-slots ((cost2 cost) (matched2 matched) (delta2 delta) (distance2 distance)) obj2
      (let ((ext-cost1 (- (+ cost1 distance1 (abs delta1)) matched1))
	    (ext-cost2 (- (+ cost2 distance2 (abs delta2)) matched2)))
	(< ext-cost1 ext-cost2)))))

(defun gap-match-runner (pat pat-ptr txt txt-ptr &optional gap-type)
  "Given two strings each with a comparison pointer, return the number
  of consecutive matching characters. if gap-type is non-nil, looks
  across ins, del, sub edits."
  (loop
    with pat-delta = (if (member gap-type '(:sub :del)) 1 0)
    with txt-delta = (if (member gap-type '(:sub :ins)) 1 0)
    for count from 0 by 1
    for pp from (+ pat-ptr pat-delta) to (1- (length pat))
    for tp from (+ txt-ptr txt-delta) to (1- (length txt))
    while (char= (char pat pp) (char txt tp))
    finally
       (return count)))

(defun last-char (string)
  "return the last character of a string"
  (char string (1- (length string))))

(defun gap-char-run (pat pat-pos &optional (gap-char #\N))
  "Return a pointer to the last character of a run of gap-chars, or
  nil if not pat-pos does not point to a gap-char."
  (if (char= gap-char (char pat pat-pos))
      (1- (or (position-if-not (lambda (c) (char= c gap-char)) pat :start pat-pos) (length pat)))
      nil))

;; create a new partial extended by an insertion
(defmethod extend-insert ((part partial) txt)
  (with-slots (moves cost pat-pos txt-pos max-pat max-txt matched distance delta) part
    (make-instance 'partial
		   :moves (cons (vector #\i (char txt txt-pos)) (copy-seq moves)) 
		   :cost (1+ cost)
		   :pat-pos pat-pos
		   :txt-pos (1+ txt-pos)
		   :max-pat max-pat
		   :max-txt max-txt
		   :matched matched
		   :distance (1- distance)
		   :delta (1+ delta))))

;; create a new partial extended by a deletion
(defmethod extend-delete ((part partial) pat)
  (with-slots (moves cost pat-pos txt-pos max-pat max-txt matched distance delta) part
    (make-instance 'partial
		   :moves (cons (vector #\d (char pat pat-pos )) (copy-seq moves)) 
		   :cost (1+ cost)
		   :pat-pos (1+ pat-pos)
		   :txt-pos txt-pos
		   :max-pat max-pat
		   :max-txt max-txt
		   :matched matched
		   :distance distance
		   :delta (1- delta))))

;; create a new partial extended by an insertion
(defmethod extend-insert-tail ((part partial) txt)
  (with-slots (moves cost pat-pos txt-pos max-pat max-txt matched distance delta) part
    (let ((tail-length (- (length txt) txt-pos)))
      (make-instance 'partial
		     :moves (cons (vector #\i (subseq txt txt-pos)) (copy-seq moves)) 
		     :cost (+ cost tail-length)
		     :pat-pos pat-pos
		     :txt-pos (+ txt-pos tail-length)
		     :max-pat max-pat
		     :max-txt max-txt
		     :matched matched
		     :distance (- distance tail-length)
		     :delta (+ delta tail-length)))))

;; create a new partial extended by deleting the entire tail
(defmethod extend-delete-tail ((part partial) pat)
  (with-slots (moves cost pat-pos txt-pos max-pat max-txt matched distance delta) part
    (let ((tail-length (- (length pat) pat-pos)))
      (make-instance 'partial
		     :moves (cons (vector #\d (subseq pat pat-pos)) (copy-seq moves)) 
		     :cost (+ cost tail-length)
		     :pat-pos (+ pat-pos tail-length)
		     :txt-pos txt-pos
		     :max-pat max-pat
		     :max-txt max-txt
		     :matched matched
		     :distance distance
		     :delta (- delta tail-length)))))

;; create a new partial extended by a run of matches
(defmethod extend-matches ((part partial) match-run)
  (with-slots (moves cost pat-pos txt-pos max-pat max-txt matched distance delta) part
    (make-instance 'partial
		   :moves (cons (vector #\m match-run) (copy-seq moves)) 
		   :cost cost
		   :pat-pos (+ pat-pos match-run)
		   :txt-pos (+ txt-pos match-run)
		   :max-pat max-pat
		   :max-txt max-txt
		   :matched (+ matched match-run)
		   :distance (- distance match-run)
		   :delta delta)))

;; extend by a deletion, including any subsequent match run
(defmethod extend-delete-matches ((part partial) pat txt)
  (with-slots (moves cost pat-pos txt-pos max-pat max-txt matched distance delta) part
    (let ((match-run (gap-match-runner pat pat-pos txt txt-pos :del))
	  (new-partial (extend-delete part pat)))
      (if (> match-run 0)
	  (extend-matches new-partial match-run)
	  new-partial))))

;; extend by an insertion, including any subsequent match run
(defmethod extend-insert-matches ((part partial) pat txt)
  (with-slots (moves cost pat-pos txt-pos max-pat max-txt matched distance delta) part
    (let ((match-run (gap-match-runner pat pat-pos txt txt-pos :ins))
	  (new-partial (extend-insert part txt)))
      (if (> match-run 0)
	  (extend-matches new-partial match-run)
	  new-partial))))

;; extend by a substitution, including any subsequent match run
(defmethod extend-substitute-matches ((part partial) pat txt)
  (with-slots (moves cost pat-pos txt-pos max-pat max-txt matched distance delta) part
    (let ((match-run (gap-match-runner pat pat-pos txt txt-pos :sub))
	  (new-partial (make-instance 'partial
				      :moves (cons (vector #\s (char txt txt-pos)) (copy-seq moves)) 
				      :cost (1+ cost)
				      :pat-pos (1+ pat-pos)
				      :txt-pos (1+ txt-pos)
				      :max-pat max-pat
				      :max-txt max-txt
				      :matched matched
				      :distance (1- distance)
				      :delta delta)))
      (if (> match-run 0)
	  (extend-matches new-partial match-run)
	  new-partial))))

;; extend by a pattern gap
(defmethod extend-pat-gap ((part partial) txt last-gap-ptr)
  (with-slots (moves cost pat-pos txt-pos max-pat max-txt matched distance delta) part
    (let ((gap-length (1+ (- last-gap-ptr pat-pos))))
      (make-instance 'partial
		     :moves (cons (vector #\g (subseq txt txt-pos (+ txt-pos gap-length)))
				  (copy-seq moves)) 
		     :cost (+ cost gap-length)
		     :pat-pos (+ pat-pos gap-length)
		     :txt-pos (+ txt-pos gap-length)
		     :max-pat max-pat
		     :max-txt max-txt
		     :matched matched
		     :distance (- distance gap-length)
		     :delta delta))))

;; extend by a text gap
(defmethod extend-txt-gap ((part partial) pat last-gap-ptr)
  (with-slots (moves cost pat-pos txt-pos max-pat max-txt matched distance delta) part
    (let ((gap-length (1+ (- last-gap-ptr txt-pos))))
      (make-instance 'partial
		     :moves (cons (vector #\g (subseq pat pat-pos (+ pat-pos gap-length)))
				  (copy-seq moves)) 
		     :cost (+ cost gap-length)
		     :pat-pos (+ pat-pos gap-length)
		     :txt-pos (+ txt-pos gap-length)
		     :max-pat max-pat
		     :max-txt max-txt
		     :matched matched
		     :distance (- distance gap-length)
		     :delta delta))))

(defmethod extend ((part partial) pat txt)
  (with-slots (moves cost matched txt-pos max-txt pat-pos max-pat distance delta) part
    (cond ((> txt-pos max-txt) ; at the end of the txt, only deletions
	   (list (extend-delete-tail part pat)))
	  ((> pat-pos max-pat) ; at the end of the pattern, only insertions
	   (list (extend-insert-tail part txt)))
	  (t (let* ((last-pat-gap-char (gap-char-run pat pat-pos))
		    (last-txt-gap-char (gap-char-run txt txt-pos))
		    (match-run (gap-match-runner pat pat-pos txt txt-pos))
		    (match-or-sub (if (> match-run 0)
				      (extend-matches part match-run)
				      (extend-substitute-matches part pat txt)))
		    (insert (extend-insert-matches part pat txt))
		    (delete (extend-delete-matches part pat txt)))
	       (cond (last-pat-gap-char ;; always take gap-char runs
		      (list (extend-pat-gap part txt last-pat-gap-char)))
		     (last-txt-gap-char
		      (list (extend-txt-gap part pat last-txt-gap-char)))
		     (t (list delete insert match-or-sub))))))))
					

(defun trim-identical-tail (string last-n)
  "trim the last-n characters from string if they are identical"
  (let* ((trim-point (- (length string) last-n))
	 (first-char (char string trim-point)))
    (if (not (find-if (lambda (c) (char-not-equal c first-char)) string :start trim-point))
	(subseq string 0 trim-point)
	string)))

(defun trim-longest-tail (pat txt)
  "when length doesn't match, remove tail of identical chars on the longer.  return both."
  (let* ((pat-length (length pat))
	 (txt-length (length txt))
	 (length-diff (abs (- pat-length txt-length))))
    (when (or (= pat-length txt-length) (< length-diff 2))
      (return-from trim-longest-tail (values pat txt)))
    (if (< pat-length txt-length)
	(values pat (trim-identical-tail txt length-diff))
	(values (trim-identical-tail pat length-diff) txt))))

(defun align (pat txt &optional (match-length 10))
  "Walk alignment of the first match-length character of pat in text
  until a match run of match-length is seen."
  (let* ((pat-length (length pat))
	 (final-match-length (if (< pat-length match-length) pat-length match-length))
	 (pat-prefix (subseq pat 0 final-match-length)))
    (search pat-prefix txt)))

(defun make-initial-partial (pat txt)
  "Given a pattern and reference text, create an initial partial solution"
  (make-instance 'partial :moves '()
			  :cost 0
			  :matched 0
			  :pat-pos 0
			  :max-pat (1- (length pat))
			  :txt-pos 0
			  :max-txt (1- (length txt))
			  :distance (length txt)
			  :delta (- (length pat) (length txt))))

(defun shift-partial (partial shift-distance longer-string ins-or-del)
  "adjust the initial partial if the pat and text are offset (:ins or :del)"
  (let ((excess (subseq longer-string 0 shift-distance)))
    (cond  ((eq ins-or-del :ins)
	    (setf (moves partial) (list (vector #\i excess)))
	    (setf (cost partial) shift-distance)
	    (setf (txt-pos partial) shift-distance)
	    (setf (distance partial) (- (distance partial) shift-distance))
	    (setf (delta partial) (+ (delta partial) shift-distance))
	    partial)
	   ((eq ins-or-del :del)
	    (setf (moves partial) '((vector #\d excess)))
	    (setf (cost partial) shift-distance)
	    (setf (pat-pos partial) shift-distance)
	    (setf (delta partial) (- (delta partial) shift-distance))
	    partial)
	   (t (error "~S is not a valid argument (:ins or :del)" ins-or-del)))))

(defun initialize-search (pat txt)
  "trim sequences and find initial alignment, return trimmed-pat
  trimmed-txt initial-partial"
  (multiple-value-bind (trimmed-pat trimmed-txt) (trim-longest-tail pat txt)
    (let ((pat-align-point (align trimmed-pat trimmed-txt))
	  (txt-align-point (align trimmed-txt trimmed-pat))
	  (initial-partial (make-initial-partial trimmed-pat trimmed-txt)))
      (unless (or pat-align-point txt-align-point)
	(return-from initialize-search (values trimmed-pat trimmed-txt initial-partial)))
      (if pat-align-point
	  (values trimmed-pat trimmed-txt (shift-partial initial-partial pat-align-point trimmed-txt :ins))
	  (values trimmed-pat trimmed-txt (shift-partial initial-partial txt-align-point trimmed-pat :del))))))

(defun lev (pat txt &optional debug)
  "compute the levenschtein or edit distance between two strings"
  (multiple-value-bind (trimmed-pat trimmed-txt initial-partial) (initialize-search pat txt)
    (let ((pq (make-instance 'heap :comparison #'compare)))
      ;; initialize the heap with an empty partial
      (shove pq initial-partial)
      (do ((best (yank pq) (yank pq))
	   (i 0 (+ i 1)))
	  ((zerop (+ (distance best) (abs (delta best)))) best)
	(when debug (break))
	(format t "iteration: ~s~%best: ~s~%" i best)
	(dolist (partial (extend best trimmed-pat trimmed-txt) nil)
	  (shove pq partial))))))

(defun rfasta (fasta-file)
  "read a fasta file and returns the header and sequence as separate string values"
  (with-open-file (fasta-stream fasta-file :direction :input)
    (let ((header (read-line fasta-stream))
	  (sequence ""))
      (loop for line = (read-line fasta-stream nil 'eof)
	    until (eq line 'eof)
	    do (setf sequence (concatenate 'string sequence line)))
      (values header sequence))))

;; read in some sequences for test data
(defvar ma-h)
(defvar ma-seq)
(multiple-value-setq (ma-h ma-seq) (rfasta "./test_data/nc2mass.fasta"))
(defvar ca-h)
(defvar ca-seq)
(multiple-value-setq (ca-h ca-seq) (rfasta "./test_data/ncov2ca.fasta"))
(defvar ref-h)
(defvar ref-seq)
(multiple-value-setq (ref-h ref-seq) (rfasta "./test_data/ncov2ref.fasta"))
