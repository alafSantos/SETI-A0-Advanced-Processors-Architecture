      .data
A:    .word 2
B:    .word 3
C:    .word 0

      .text
main:
      ld r4,A(r0)
      ld r5,B(r0)
      dmult r4,r5
      mflo r3
      sd r3, C(r0)
      halt          
