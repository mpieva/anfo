; Everything in here was copied nearly verbatim from
; http://srfi.schemers.org/srfi-88/srfi-88.html, therefore:
;
; Copyright (C) Marc Feeley (2006). All Rights Reserved. 
; Permission is hereby granted, free of charge, to any person obtaining
; a copy of this software and associated documentation files (the
; "Software"), to deal in the Software without restriction, including
; without limitation the rights to use, copy, modify, merge, publish,
; distribute, sublicense, and/or sell copies of the Software, and to
; permit persons to whom the Software is furnished to do so, subject to
; the following conditions: 
;
; The above copyright notice and this permission notice shall be
; included in all copies or substantial portions of the Software. 
;
; THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
; EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
; MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
; IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
; CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
; TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
; SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

(define real-symbol? symbol?)
    (define real-symbol->string symbol->string)
    (define real-string->symbol string->symbol)

    (define looks-like-an-unquoted-keyword?
      (lambda (s)
        (let ((n (string-length s)))
          (and (> n 1)
               (char=? (string-ref s (- n 1)) #\:)))))

    (set! symbol?
      (lambda (obj)
        (and (real-symbol? obj)
             (not (looks-like-an-unquoted-keyword?
                   (real-symbol->string obj))))))

    (define keyword?
      (lambda (obj)
        (and (real-symbol? obj)
             (looks-like-an-unquoted-keyword?
              (real-symbol->string obj)))))

    (set! symbol->string real-symbol->string)

    (define keyword->string
      (lambda (k)
        (let* ((s (real-symbol->string k))
               (n (string-length s)))
          (substring s 0 (- n 1))))) ; remove the colon

    (set! string->symbol
      (lambda (s)
        (if (looks-like-an-unquoted-keyword? s)
            (error "sorry... the symbol would look like a keyword!")
            (real-string->symbol s))))

    (define string->keyword
      (lambda (s)
        (let ((s-colon (string-append s ":")))
          (if (looks-like-an-unquoted-keyword? s-colon)
              (real-string->symbol s-colon)
              (error "sorry... the keyword would look like a symbol!")))))

(provide 'srfi-88)
