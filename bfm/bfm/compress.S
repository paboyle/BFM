.text
.align 2
.globl vmx_compress
.type vmx_compress,@function
vmx_compress:
		la   %r1,	-960(%r1)
	std   %r0,	 256(%r1)
	std   %r2,	 272(%r1)
	std   %r3,	 280(%r1)
	std   %r4,	 288(%r1)
	std   %r5,	 296(%r1)
	std   %r6,	 304(%r1)
	std   %r7,	 312(%r1)
	std   %r8,	 320(%r1)
	std   %r9,	 328(%r1)
	std   %r10,	 336(%r1)
	std   %r11,	 344(%r1)
	std   %r12,	 352(%r1)
	std   %r13,	 360(%r1)
	std   %r14,	 368(%r1)
	std   %r15,	 376(%r1)
	std   %r16,	 384(%r1)
	std   %r17,	 392(%r1)
	std   %r18,	 400(%r1)
	std   %r19,	 408(%r1)
	std   %r20,	 416(%r1)
	std   %r21,	 424(%r1)
	std   %r22,	 432(%r1)
	std   %r23,	 440(%r1)
	std   %r24,	 448(%r1)
	std   %r25,	 456(%r1)
	std   %r26,	 464(%r1)
	std   %r27,	 472(%r1)
	std   %r28,	 480(%r1)
	std   %r29,	 488(%r1)
	std   %r30,	 496(%r1)
	std   %r31,	 504(%r1)
	stfd  0,	   0(%r1)
	stfd  1,	   8(%r1)
	stfd  2,	  16(%r1)
	stfd  3,	  24(%r1)
	stfd  4,	  32(%r1)
	stfd  5,	  40(%r1)
	stfd  6,	  48(%r1)
	stfd  7,	  56(%r1)
	stfd  8,	  64(%r1)
	stfd  9,	  72(%r1)
	stfd  10,	  80(%r1)
	stfd  11,	  88(%r1)
	stfd  12,	  96(%r1)
	stfd  13,	 104(%r1)
	stfd  14,	 112(%r1)
	stfd  15,	 120(%r1)
	stfd  16,	 128(%r1)
	stfd  17,	 136(%r1)
	stfd  18,	 144(%r1)
	stfd  19,	 152(%r1)
	stfd  20,	 160(%r1)
	stfd  21,	 168(%r1)
	stfd  22,	 176(%r1)
	stfd  23,	 184(%r1)
	stfd  24,	 192(%r1)
	stfd  25,	 200(%r1)
	stfd  26,	 208(%r1)
	stfd  27,	 216(%r1)
	stfd  28,	 224(%r1)
	stfd  29,	 232(%r1)
	stfd  30,	 240(%r1)
	stfd  31,	 248(%r1)
	li   %r18,   0
	li   %r19,  32
	li   %r20,  64
	li   %r21,  96
	li   %r22, 128
	li   %r23, 160
	or    %r11, %r3, %r3 /*number of half spinors*/
	or    %r12, %r4, %r4 /*input*/
	or    %r13, %r5, %r5 /*output*/
	mtctr 	%r11
_vmx_compress_lab1:
	/*Load a half spinor (3x2) complex simd*/
	qvlfdx   0,%r18,%r12
	qvlfdx   1,%r19,%r12
	qvlfdx   2,%r20,%r12
	qvlfdx   3,%r21,%r12
	qvlfdx   4,%r22,%r12
	qvlfdx   5,%r23,%r12
	/*Find the absolute value*/
	qvfabs   6,0
	qvfabs   7,1
	qvfabs   8,2
	qvfabs   9,3
	qvfabs   10,4
	qvfabs   11,5
	/*Find the absolute max*/
	qvfcmpgt 12,6,7
	qvfcmpgt 13,8,9
	qvfcmpgt 14,10,11
	qvfsel   6,12,6,7
	qvfsel   7,13,8,9
	qvfsel   8,14,10,11
	qvfcmpgt 12,6,7
	qvfsel   6,12,6,7
	qvfcmpgt 12,6,8
	qvfsel   6,12,6,8
	qvfaligni 7,6,6,1
	qvfcmpgt 12,6,7
	qvfsel   6,12,6,7
	qvfaligni 7,6,6,2
	qvfcmpgt 12,6,7
	qvfsel   6,12,6,7
	qvfaligni 7,6,6,3
	qvfcmpgt 12,6,7
	qvfsel   6,12,6,7
	/*Divide by the max use qvfre, qvfmul*/

	/*multiply by 64 =2^6*/
	
	/*convert to integer*/

	/*store need to get to bytes*/
	
	qvstfdx  6,%r18,%r12
	qvstfdx  6,%r19,%r12
	qvstfdx  6,%r20,%r12
	qvstfdx  6,%r21,%r12
	qvstfdx  6,%r22,%r12
	qvstfdx  6,%r23,%r12
	bdnz _vmx_compress_lab1
	lfd   31,	 248(%r1)
	lfd   30,	 240(%r1)
	lfd   29,	 232(%r1)
	lfd   28,	 224(%r1)
	lfd   27,	 216(%r1)
	lfd   26,	 208(%r1)
	lfd   25,	 200(%r1)
	lfd   24,	 192(%r1)
	lfd   23,	 184(%r1)
	lfd   22,	 176(%r1)
	lfd   21,	 168(%r1)
	lfd   20,	 160(%r1)
	lfd   19,	 152(%r1)
	lfd   18,	 144(%r1)
	lfd   17,	 136(%r1)
	lfd   16,	 128(%r1)
	lfd   15,	 120(%r1)
	lfd   14,	 112(%r1)
	lfd   13,	 104(%r1)
	lfd   12,	  96(%r1)
	lfd   11,	  88(%r1)
	lfd   10,	  80(%r1)
	lfd   9,	  72(%r1)
	lfd   8,	  64(%r1)
	lfd   7,	  56(%r1)
	lfd   6,	  48(%r1)
	lfd   5,	  40(%r1)
	lfd   4,	  32(%r1)
	lfd   3,	  24(%r1)
	lfd   2,	  16(%r1)
	lfd   1,	   8(%r1)
	lfd   0,	   0(%r1)
	ld   %r31,	 504(%r1)
	ld   %r30,	 496(%r1)
	ld   %r29,	 488(%r1)
	ld   %r28,	 480(%r1)
	ld   %r27,	 472(%r1)
	ld   %r26,	 464(%r1)
	ld   %r25,	 456(%r1)
	ld   %r24,	 448(%r1)
	ld   %r23,	 440(%r1)
	ld   %r22,	 432(%r1)
	ld   %r21,	 424(%r1)
	ld   %r20,	 416(%r1)
	ld   %r19,	 408(%r1)
	ld   %r18,	 400(%r1)
	ld   %r17,	 392(%r1)
	ld   %r16,	 384(%r1)
	ld   %r15,	 376(%r1)
	ld   %r14,	 368(%r1)
	ld   %r13,	 360(%r1)
	ld   %r12,	 352(%r1)
	ld   %r11,	 344(%r1)
	ld   %r10,	 336(%r1)
	ld   %r9,	 328(%r1)
	ld   %r8,	 320(%r1)
	ld   %r7,	 312(%r1)
	ld   %r6,	 304(%r1)
	ld   %r5,	 296(%r1)
	ld   %r4,	 288(%r1)
	ld   %r3,	 280(%r1)
	ld   %r2,	 272(%r1)
	ld   %r0,	 256(%r1)
	la   %r1,	 960(%r1)
	blr
.globl _vmx_compress_end
_vmx_compress_end:
.size vmx_compress, .-vmx_compress
	
