	;; B1=abs(A1) 
	.data
A1:	.double -10.0
B1:	.double 0.0
zero:	.double 0.0
	.text
main:
	
	l.d f0,zero(r0)
	l.d f2, A1(r0)	        ; f2=A1
	mov.d f4, f2		; f4=B1 = A1
	c.lt.d 7, f2, f0	; B1 < 0 ?
	bc1f 7, done		; if false, nothing to do
	sub.d f4, f0, f2	; f4=-A=|A|= 0.0 - f4
done:	
	s.d f4,B1(r0)
	halt          
 
