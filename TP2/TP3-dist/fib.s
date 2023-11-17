; Compute first twelve Fibonacci numbers and put in array, then print
      .data
fibs: .space 48          ; "array" of 12 32_bits words to contain fib values
size: .word  12          ; size of "array" 
      .text
      addi $t0, r0, fibs    ; load address of array
      lw   $t5, size(r0)    ; load array size
      addi $t2, r0,1        ; 1 is first and second Fib. number
      sw   $t2, 0($t0)      ; F[0] = 1
      sw   $t2, 4($t0)      ; F[1] = F[0] = 1
      addi $t1, $t5, -2     ; Counter for loop, will execute (size-2) times
loop: lw   $t3, 0($t0)      ; Get value from array F[n] -> $t3
      lw   $t4, 4($t0)      ; Get value from array F[n+1] -> *t4
      add  $t2, $t3, $t4    ; $t2 = F[n] + F[n+1]
      sw   $t2, 8($t0)      ; Store F[n+2] = F[n] + F[n+1] in array
      addi $t0, $t0, 4      ; increment address of Fib. number source
      addi $t1, $t1, -1     ; decrement loop counter
      bnez $t1, loop        ; repeat if not finished yet.
      nop
      halt
	
     
