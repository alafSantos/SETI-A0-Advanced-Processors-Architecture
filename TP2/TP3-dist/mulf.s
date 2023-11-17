      .data
A:    .double 10.0
B:    .double 8.0
C:    .double 0.0

      .text
main:
      l.d f0,A(r0)
      l.d f2,B(r0)
      mul.d f4,f0,f2
      s.d f4,C(r0)
      halt          
