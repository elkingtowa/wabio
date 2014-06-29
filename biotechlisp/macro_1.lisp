(defmacro nil! (x)
  (list 'setf x nil))
