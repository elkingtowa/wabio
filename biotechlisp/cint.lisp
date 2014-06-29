#|
Date: Mon, 10 Dec 2001 18:16:48 -0500 (EST)
From: Guy Steele
To: ll1-discuss@ai.mit.edu

This is an interpreter for a tiny Scheme-like language
(I'll call it FOO) with just the following constructs:

        numeric literals
        variable names
        LAMBDA expressions
        IF
        function calls

and three built-in functions:

        +
        *
        CALL/CC

(I took the cheesy way out on defining the functions:
I just made undefined variables evaluate to themselves
and have the @apply function check for those names.
That's okay because the data domain of this language
is just numbers---there aren't any operations on symbols.)

I spelled the name of each function with a leading "@" purely
to avoid conflict with the Common Lisp functions of the
same name.

You can evaluate an expression by typing at the Common Lisp
top level:

        (@eval '<expression> '() #'(lambda (x) x))

Note two things about this piece of code:

(1) Every call from one Common Lisp function to another is
a tail-call.  In other words, in effect I am not using the
Common Lisp stack at all to keep information about the state
of the FOO program being interpreted.

(2) Every LAMBDA expression in the code of the interpreter
is a continuation: it says what to do next when the call to
any given @-routine is "finished".

(3) Every @-routine takes a continuation "cont" and always
finishes either by calling cont (as a tail-call) or by
calling another @-routine (as a tail-call).

If you keep in mind that "#'" means roughly "allocate  (in the heap)
a closure for the following LAMBDA expression" and that a closure can
refer to lexical variables visible to the LAMBDA expression, you can see
that one continuation can know about another, which knows about another,
and so on; this chain is sometimes called the "control stack", but in
this implementation it's all in the heap.

------------------------------------------------------------------

Given this (pathological) test program

((call/cc
  (lambda (goto)
    (letrec ((start
              (lambda ()
                (print "start")
                (goto next)))
             (froz
              (lambda ()
                (print "froz")
                (goto last)))
             (next
              (lambda ()
                (print "next")
                (goto froz)))
             (last
              (lambda ()
                (print "last")
                (+ 3 4))))
      start))))

the output is:

start
next
froz
last
7
|#



(defun @eval (exp env cont)
  (cond ((numberp exp) (funcall cont exp))
        ((stringp exp) (funcall cont exp))
        ((symbolp exp) (@lookup exp env cont))
        ((eq (first exp) 'LAMBDA)
         (funcall cont (list 'CLOSURE (second exp) (rest (rest exp)) env)))
        ((eq (first exp) 'IF)
         (@eval (second exp) env
                #'(lambda (test)
                    (@eval (cond (test (second exp)) (t (third exp))) env cont))))
        ((eq (first exp) 'LETREC)
         (let ((newenv (pairlis (mapcar #'first (second exp))
                                (make-list (length (second exp)))
                                env)))
           (@evletrec (second exp) newenv (third exp) newenv cont)))
        (t (@eval (first exp) env
                  #'(lambda (fn)
                      (@evlis (rest exp) env
                              #'(lambda (args) (@apply fn args cont))))))))

(defun @lookup (name env cont)
  (cond ((null env) (funcall cont name))
        ((eq (car (first env)) name) (funcall cont (cdr (first env))))
        (t (@lookup name (rest env) cont))))

(defun @evlis (exps env cont)
  (cond ((null exps) (funcall cont '()))
        (t (@eval (first exps) env
                  #'(lambda (arg)
                      (@evlis (rest exps) env
                              #'(lambda (args) (funcall cont (cons arg args)))))))))

(defun @evletrec (bindings slots body env cont)
  (cond ((null bindings) (@eval body env cont))
        (t (@eval (second (first bindings)) env
                  #'(lambda (fn)
                      (rplacd (first slots) fn)  ;the side effect that "ties the knot"
                      (@evletrec (rest bindings) (rest slots) body env cont))))))

(defun @apply (fn args cont)
  (cond ((eq fn '+) (funcall cont (+ (first args) (second args))))
        ((eq fn '*) (funcall cont (* (first args) (second args))))
        ((eq fn 'print)
         (princ (first args))
         (fresh-line)
         (funcall cont (first args)))
        ((eq fn 'call/cc)
         (@apply (first args) (list (list 'CONTINUATION cont)) cont))
        ((atom fn) (funcall cont 'UNDEFINED-FUNCTION))
        ((eq (first fn) 'CLOSURE)
         (@evlis (third fn) (pairlis (second fn) args (fourth fn))
                 #'(lambda (vals) (funcall cont (first (last vals))))))
        ((eq (first fn) 'CONTINUATION)
         (funcall (second fn) (first args)))
        (t (funcall cont 'UNDEFINED-FUNCTION))))


