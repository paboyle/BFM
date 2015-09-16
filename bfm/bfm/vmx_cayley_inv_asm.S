/* 
 * BAGEL machine generated output.  
 *   
 * It is provided under the GNU pubic License V2  
 * It is provided as is, and is not guaranteed fit for any purpose.
 * BAGEL was written by Peter Boyle  
 */  
.text
.align 16,0x90
.globl vmx_cayley_inv
vmx_cayley_inv:
	pushq    %r10
	movq     %rsp, %r10
	andq     $-64, %rsp
	movq     %r10, (%rsp)
	movq     $0xf000cc000100aa, %r10
	kextract $0, %r10, %k6   # even elements
	kextract $2, %r10, %k4   # even pairs
	kextract $1, %r10, %k5   # scalar flop
	kextract $3, %r10, %k2   # scalar flop
	knot     %k6,  %k7       # odd elements
	knot     %k4,  %k3       # odd pairs
	# ----- end of Enter_Routine -----

	addq          $-2176, %rsp
	movq          %rdi, %r11
	movq          %rbx, 2112(%rsp)
	movq          %rbp, 2120(%rsp)
	movq          %r12, 2144(%rsp)
	movq          %r13, 2152(%rsp)
	movq          %r14, 2160(%rsp)
	movq          %r15, 2168(%rsp)
	movq             0(%r11), %r8
	addq          $8, %r11
	movq             0(%r11), %rcx
	addq          $8, %r11
	movq             0(%r11), %r9
	addq          $8, %r11
	movq             0(%r11), %rbp
	addq          $8, %r11
	movq             0(%r11), %rax
	addq          $8, %r11
_vmx_cayley_inv_lab0:
	movq          %r8, %r12
	movq          %rcx, %r13
	movq          %rax, %rbx
	vmovapd          0(%r12), %zmm12
	vbroadcastsd     0(%rbx), %zmm24
	vbroadcastsd     8(%rbx), %zmm25
	vmovapd        384(%r12), %zmm18
	vmovapd        128(%r12), %zmm14
	vmovapd        512(%r12), %zmm20
	vmulpd        %zmm18, %zmm25, %zmm6
	vmovapd         64(%r12), %zmm13
	vmovapd        448(%r12), %zmm19
	vmulpd        %zmm20, %zmm25, %zmm8
	vmovapd        192(%r12), %zmm15
	vmovapd        576(%r12), %zmm21
	vmulpd        %zmm19, %zmm25, %zmm7
	vmovapd        320(%r12), %zmm17
	vmovapd        704(%r12), %zmm23
	vmulpd        %zmm21, %zmm25, %zmm9
	vmovapd        256(%r12), %zmm16
	vmovapd        640(%r12), %zmm22
	vmulpd        %zmm23, %zmm25, %zmm11
	vmovnrngoapd  %zmm12,    0(%r13)
	movq          %rbp, %r10
	vmulpd        %zmm22, %zmm25, %zmm10
	vmovnrngoapd  %zmm13,   64(%r13)
	addq          $-1, %r10
	vmovnrngoapd  %zmm14,  128(%r13)
	vmovnrngoapd  %zmm15,  192(%r13)
	vmovnrngoapd  %zmm16,  256(%r13)
	vmovnrngoapd  %zmm17,  320(%r13)
	vmovnrngoapd  %zmm18,  384(%r13)
	vmovnrngoapd  %zmm19,  448(%r13)
	vmovnrngoapd  %zmm20,  512(%r13)
	vmovnrngoapd  %zmm21,  576(%r13)
	vmovnrngoapd  %zmm22,  640(%r13)
	vmovnrngoapd  %zmm23,  704(%r13)
_vmx_cayley_inv_lab1:
	addq          $768, %r12
	addq          $768, %r13
	addq          $16, %rbx
	vmovapd        384(%r12), %zmm18
	vbroadcastsd     0(%rbx), %zmm24
	vbroadcastsd     8(%rbx), %zmm25
	vmovapd          0(%r12), %zmm0
	vmovapd        576(%r12), %zmm21
	vfmadd231pd   %zmm18, %zmm25, %zmm6
	vmovapd        192(%r12), %zmm3
	vfmadd213pd   %zmm0, %zmm24, %zmm12
	vmovapd        448(%r12), %zmm19
	vfmadd231pd   %zmm21, %zmm25, %zmm9
	vmovapd         64(%r12), %zmm1
	vfmadd213pd   %zmm3, %zmm24, %zmm15
	vmovapd        640(%r12), %zmm22
	vfmadd231pd   %zmm19, %zmm25, %zmm7
	vmovapd        256(%r12), %zmm4
	vfmadd213pd   %zmm1, %zmm24, %zmm13
	vmovapd        512(%r12), %zmm20
	vfmadd231pd   %zmm22, %zmm25, %zmm10
	vmovapd        128(%r12), %zmm2
	vfmadd213pd   %zmm4, %zmm24, %zmm16
	vmovapd        704(%r12), %zmm23
	vfmadd231pd   %zmm20, %zmm25, %zmm8
	vmovapd        320(%r12), %zmm5
	vfmadd213pd   %zmm2, %zmm24, %zmm14
	vmovnrngoapd  %zmm12,    0(%r13)
	vfmadd231pd   %zmm23, %zmm25, %zmm11
	addq          $-1, %r10
	vfmadd213pd   %zmm5, %zmm24, %zmm17
	vmovnrngoapd  %zmm13,   64(%r13)
	orq           %r10, %r10
	vmovnrngoapd  %zmm14,  128(%r13)
	vmovnrngoapd  %zmm15,  192(%r13)
	vmovnrngoapd  %zmm16,  256(%r13)
	vmovnrngoapd  %zmm17,  320(%r13)
	vmovnrngoapd  %zmm18,  384(%r13)
	vmovnrngoapd  %zmm19,  448(%r13)
	vmovnrngoapd  %zmm20,  512(%r13)
	vmovnrngoapd  %zmm21,  576(%r13)
	vmovnrngoapd  %zmm22,  640(%r13)
	vmovnrngoapd  %zmm23,  704(%r13)
	jg     _vmx_cayley_inv_lab1
	addq          $16, %rbx
	movq          %rbp, %r10
	vbroadcastsd     0(%rbx), %zmm26
	addq          $8, %rbx
	addq          $-1, %r10
	vmulpd        %zmm6, %zmm26, %zmm6
	vmulpd        %zmm7, %zmm26, %zmm7
	vmulpd        %zmm8, %zmm26, %zmm8
	vmulpd        %zmm9, %zmm26, %zmm9
	vmulpd        %zmm10, %zmm26, %zmm10
	vmovnrngoapd  %zmm6,  384(%r13)
	vmulpd        %zmm11, %zmm26, %zmm11
	vmulpd        %zmm12, %zmm26, %zmm0
	vmovnrngoapd  %zmm7,  448(%r13)
	vmulpd        %zmm13, %zmm26, %zmm1
	vmulpd        %zmm14, %zmm26, %zmm2
	vmovnrngoapd  %zmm8,  512(%r13)
	vmulpd        %zmm15, %zmm26, %zmm3
	vmulpd        %zmm16, %zmm26, %zmm4
	vmovnrngoapd  %zmm0,    0(%r13)
	vmulpd        %zmm17, %zmm26, %zmm5
	vmovnrngoapd  %zmm1,   64(%r13)
	vmovnrngoapd  %zmm2,  128(%r13)
	vmovnrngoapd  %zmm3,  192(%r13)
	vmovnrngoapd  %zmm4,  256(%r13)
	vmovnrngoapd  %zmm5,  320(%r13)
	vmovnrngoapd  %zmm9,  576(%r13)
	vmovnrngoapd  %zmm10,  640(%r13)
	vmovnrngoapd  %zmm11,  704(%r13)
_vmx_cayley_inv_lab2:
	vbroadcastsd     0(%rbx), %zmm26
	vbroadcastsd     8(%rbx), %zmm24
	vbroadcastsd    16(%rbx), %zmm25
	addq          $24, %rbx
	addq          $-768, %r13
	addq          $-1, %r10
	vmovapd          0(%r13), %zmm12
	vmovapd        384(%r13), %zmm18
	vmovapd        192(%r13), %zmm15
	vmulpd        %zmm12, %zmm26, %zmm12
	vmovapd        576(%r13), %zmm21
	vmulpd        %zmm18, %zmm26, %zmm18
	vmovapd         64(%r13), %zmm13
	vmulpd        %zmm15, %zmm26, %zmm15
	vmovapd        448(%r13), %zmm19
	vmulpd        %zmm21, %zmm26, %zmm21
	vmovapd        256(%r13), %zmm16
	vmulpd        %zmm13, %zmm26, %zmm13
	vmovapd        640(%r13), %zmm22
	vmulpd        %zmm19, %zmm26, %zmm19
	vmovapd        128(%r13), %zmm14
	vmulpd        %zmm16, %zmm26, %zmm16
	vmovapd        512(%r13), %zmm20
	vmulpd        %zmm22, %zmm26, %zmm22
	vmovapd        320(%r13), %zmm17
	vmulpd        %zmm14, %zmm26, %zmm14
	vmovapd        704(%r13), %zmm23
	vmulpd        %zmm20, %zmm26, %zmm20
	orq           %r10, %r10
	vmulpd        %zmm17, %zmm26, %zmm17
	vmulpd        %zmm23, %zmm26, %zmm23
	vfmadd231pd   %zmm0, %zmm25, %zmm12
	vfmadd231pd   %zmm3, %zmm25, %zmm15
	vfmadd231pd   %zmm1, %zmm25, %zmm13
	vfmadd231pd   %zmm4, %zmm25, %zmm16
	vfmadd231pd   %zmm2, %zmm25, %zmm14
	vmovnrngoapd  %zmm12,    0(%r13)
	vfmadd231pd   %zmm5, %zmm25, %zmm17
	vfmadd213pd   %zmm18, %zmm24, %zmm6
	vmovnrngoapd  %zmm13,   64(%r13)
	vfmadd213pd   %zmm21, %zmm24, %zmm9
	vfmadd213pd   %zmm19, %zmm24, %zmm7
	vmovnrngoapd  %zmm14,  128(%r13)
	vfmadd213pd   %zmm22, %zmm24, %zmm10
	vfmadd213pd   %zmm20, %zmm24, %zmm8
	vmovnrngoapd  %zmm6,  384(%r13)
	vfmadd213pd   %zmm23, %zmm24, %zmm11
	vmovnrngoapd  %zmm7,  448(%r13)
	vmovnrngoapd  %zmm8,  512(%r13)
	vmovnrngoapd  %zmm9,  576(%r13)
	vmovnrngoapd  %zmm10,  640(%r13)
	vmovnrngoapd  %zmm11,  704(%r13)
	vmovnrngoapd  %zmm15,  192(%r13)
	vmovnrngoapd  %zmm16,  256(%r13)
	vmovnrngoapd  %zmm17,  320(%r13)
	jg     _vmx_cayley_inv_lab2
	imulq         $768, %rbp, %r12
	addq          $-1, %r9
	addq          %r12, %r8
	addq          %r12, %rcx
	orq           %r9, %r9
	jg     _vmx_cayley_inv_lab0
	movq          2168(%rsp), %r15
	movq          2160(%rsp), %r14
	movq          2152(%rsp), %r13
	movq          2144(%rsp), %r12
	movq          2120(%rsp), %rbp
	movq          2112(%rsp), %rbx
	addq          $2176, %rsp

	# ----- Exit_Routine -----
	movq     (%rsp), %rsp
	popq     %r10
	ret
.align 16,0x90
.type vmx_cayley_inv,@function
.size vmx_cayley_inv,.-vmx_cayley_inv
.data
