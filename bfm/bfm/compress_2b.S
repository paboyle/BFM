.text
.align 2
.globl vmx_compress2b
.type vmx_compress2b,@function
vmx_compress2b:
	la   %r1,	-960(%r1)
	std   %r13,	 360(%r1)
	std   %r14,	 368(%r1)
	std   %r15,	 376(%r1)
	std   %r16,	 384(%r1)
	std   %r17,	 392(%r1)
	rotldi  %r3,%r3,1
	
	li %r17,0
	ori %r17,%r17,0xFFFF
	sldi %r17,%r17,32
	ori %r17,%r17,0xFFFF
	sldi %r17,%r17,16

	li %r6,32
	li %r7,64
	mtctr 	%r3
_vmx_compress2b_lab1:
	ld    %r10,0(%r4)
	ld    %r11,8(%r4)
	ld    %r12,16(%r4)
	ld    %r13,24(%r4)
	ld    %r14,32(%r4)
	ld    %r15,40(%r4)
	addi %r4,%r4,48

	and   %r10,%r17,%r10
	dcbt %r6,%r4
	and   %r11,%r17,%r11
	dcbt %r7,%r4
	and   %r12,%r17,%r12
	and   %r13,%r17,%r13
	and   %r14,%r17,%r14
	and   %r15,%r17,%r15

	srdi %r11,%r11,16
	srdi %r13,%r13,16
	srdi %r15,%r15,16
	
	or   %r10,%r10,%r11
	or   %r12,%r12,%r13
	std   %r10,0(%r5)
	or   %r14,%r14,%r15
	std   %r12,8(%r5)
	std   %r14,16(%r5)
	addi %r5,%r5,24
	
	bdnz _vmx_compress2b_lab1
	
	ld   %r17,	 392(%r1)
	ld   %r16,	 384(%r1)
	ld   %r15,	 376(%r1)
	ld   %r14,	 368(%r1)
	ld   %r13,	 360(%r1)
	la   %r1,	 960(%r1)
	blr
.globl _vmx_compress2b_end
_vmx_compress2b_end:
.size vmx_compress2b, .-vmx_compress2b
	
.text
.align 2
.globl vmx_decompress2b
.type vmx_decompress2b,@function
vmx_decompress2b:
	la   %r1,	-960(%r1)
	std   %r13,	 360(%r1)
	std   %r14,	 368(%r1)
	std   %r15,	 376(%r1)
	std   %r16,	 384(%r1)
	rotldi %r3,%r3,1
	mtctr 	%r3

	li %r16,0
	ori %r16,%r16,0xFFFF
	rotldi %r16,%r16,32
	ori %r16,%r16,0xFFFF
	li %r6,32

_vmx_decompress2b_lab1:

	ld %r10,0(%r4)  
	ld %r11,8(%r4)
	ld %r12,16(%r4)

	and  %r13,%r10,%r16
	andc %r10,%r10,%r16
	sldi %r13,%r13,16
	std %r10,0(%r5)
	std %r13,8(%r5)
	
	and  %r14,%r11,%r16
	andc %r11,%r11,%r16
	sldi %r14,%r14,16
	std %r11,16(%r5)
	std %r14,24(%r5)

	and  %r15,%r12,%r16
	andc %r12,%r12,%r16
	sldi %r15,%r15,16
	std %r12,32(%r5)
	std %r15,40(%r5)

	addi  %r4,%r4,24
	addi  %r5,%r5,48
	dcbt %r6,%r4
	
	bdnz _vmx_decompress2b_lab1
	
	ld   %r16,	 384(%r1)
	ld   %r15,	 376(%r1)
	ld   %r14,	 368(%r1)
	ld   %r13,	 360(%r1)
	la   %r1,	 960(%r1)
	blr
