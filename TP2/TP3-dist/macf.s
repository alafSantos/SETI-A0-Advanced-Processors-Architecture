      .data
A:    .double 10.0
B:    .double 8.0
C:    .double 15.0
D: 	.double 0.0

      .text
main:
      l.d f0,A(r0)
      l.d f1,B(r0)
	 l.d f2,C(r0)
      mul.d f1,f0,f1
	add.d f2,f2,f1
      s.d f2,D(r0)
      halt          
