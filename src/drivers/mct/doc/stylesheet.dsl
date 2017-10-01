<!DOCTYPE style-sheet PUBLIC "-//James Clark//DTD DSSSL Style Sheet//EN" [
<!ENTITY html-ss
  PUBLIC "-//Norman Walsh//DOCUMENT DocBook HTML Stylesheet//EN" CDATA dsssl>
]>

<style-specification id="html" use="html-stylesheet">
<style-specification-body>

;;Default extension for filenames
(define %html-ext% ".html")
;;What font would you like for the body?
(define %body-font-family%
 "Arial")

(element emphasis
(if (equal? (normalize "bold") (attribute-string (normalize "role")))
    ($bold-seq$)
    ($italic-seq$)))

(element tgroup
  (let* ((wrapper   (parent (current-node)))
	  (frameattr (attribute-string (normalize "frame") wrapper))
	   (pgwide    (attribute-string (normalize "pgwide") wrapper))
	    (footnotes (select-elements (descendants (current-node))
					     (normalize "footnote")))
	     (border (if (equal? frameattr (normalize "none"))
			      '(("BORDER" "0"))
			           '(("BORDER" "1"))))
	      (bgcolor '(("BGCOLOR" "#E0E0E0")))
	       (width (if (equal? pgwide "1")
			      (list (list "WIDTH" ($table-width$)))
			          '()))
	        (head (select-elements (children (current-node)) (normalize "thead")))
		 (body (select-elements (children (current-node)) (normalize "tbody")))
		  (feet (select-elements (children (current-node)) (normalize "tfoot"))))
    (make element gi: "TABLE"
	    attributes: (append
			        border
				       width
				              bgcolor
					             '(("CELLSPACING" "0"))
						            '(("CELLPADDING" "4"))
							           (if %cals-table-class%
								          (list (list "CLASS" %cals-table-class%))
									     '()))
	      (process-node-list head)
	        (process-node-list body)
		  (process-node-list feet)
		    (make-table-endnotes))))


;;Should verbatim items be 'shaded' with a table?
(define %shade-verbatim%
 #t)

;;Define shade-verbatim attributes
(define ($shade-verbatim-attr$)
 (list
  (list "BORDER" "0")
  (list "BGCOLOR" "#E0E0E0")
  (list "WIDTH" ($table-width$))))


</style-specification-body>
</style-specification>
<external-specification id="html-stylesheet" document="html-ss">
</style-sheet>

