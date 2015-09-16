/* 
 * BAGEL machine generated output.  
 *   
 * It is provided under the GNU pubic License V2  
 * It is provided as is, and is not guaranteed fit for any purpose.
 * BAGEL was written by Peter Boyle  
 */  
.text
.align 16,0x90
.globl vmx_cayley_inv_s
vmx_cayley_inv_s:
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
_vmx_cayley_inv_s_lab0:
	movq          %r8, %r12
	movq          %rcx, %r13
	movq          %rax, %rbx
	vmovaps          0(%r12), %zmm12
	vbroadcastsd     0(%rbx), %zmm24
	vbroadcastsd     8(%rbx), %zmm25
	vmovaps        192(%r12), %zmm18
	vmovaps         64(%r12), %zmm14
	vmovaps        256(%r12), %zmm20
	vmulpd        %zmm18, %zmm25, %zmm6
	vmovaps         32(%r12), %zmm13
	vmovaps        224(%r12), %zmm19
	vmulpd        %zmm20, %zmm25, %zmm8
	vmovaps         96(%r12), %zmm15
	vmovaps        288(%r12), %zmm21
	vmulpd        %zmm19, %zmm25, %zmm7
	vmovaps        160(%r12), %zmm17
	vmovaps        352(%r12), %zmm23
	vmulpd        %zmm21, %zmm25, %zmm9
	vmovaps        128(%r12), %zmm16
	vmovaps        320(%r12), %zmm22
	vmulpd        %zmm23, %zmm25, %zmm11
	vmovnrngoaps  %zmm12,    0(%r13)
	movq          %rbp, %r10
	vmulpd        %zmm22, %zmm25, %zmm10
	vmovnrngoaps  %zmm13,   32(%r13)
	addq          $-1, %r10
	vmovnrngoaps  %zmm14,   64(%r13)
	vmovnrngoaps  %zmm15,   96(%r13)
	vmovnrngoaps  %zmm16,  128(%r13)
	vmovnrngoaps  %zmm17,  160(%r13)
	vmovnrngoaps  %zmm18,  192(%r13)
	vmovnrngoaps  %zmm19,  224(%r13)
	vmovnrngoaps  %zmm20,  256(%r13)
	vmovnrngoaps  %zmm21,  288(%r13)
	vmovnrngoaps  %zmm22,  320(%r13)
	vmovnrngoaps  %zmm23,  352(%r13)
_vmx_cayley_inv_s_lab1:
	addq          $384, %r12
	addq          $384, %r13
	addq          $16, %rbx
	vmovaps        192(%r12), %zmm18
	vbroadcastsd     0(%rbx), %zmm24
	vbroadcastsd     8(%rbx), %zmm25
	vmovaps          0(%r12), %zmm0
	vmovaps        288(%r12), %zmm21
	vfmadd231pd   %zmm18, %zmm25, %zmm6
	vmovaps         96(%r12), %zmm3
	vfmadd213pd   %zmm0, %zmm24, %zmm12
	vmovaps        224(%r12), %zmm19
	vfmadd231pd   %zmm21, %zmm25, %zmm9
	vmovaps         32(%r12), %zmm1
	vfmadd213pd   %zmm3, %zmm24, %zmm15
	vmovaps        320(%r12), %zmm22
	vfmadd231pd   %zmm19, %zmm25, %zmm7
	vmovaps        128(%r12), %zmm4
	vfmadd213pd   %zmm1, %zmm24, %zmm13
	vmovaps        256(%r12), %zmm20
	vfmadd231pd   %zmm22, %zmm25, %zmm10
	vmovaps         64(%r12), %zmm2
	vfmadd213pd   %zmm4, %zmm24, %zmm16
	vmovaps        352(%r12), %zmm23
	vfmadd231pd   %zmm20, %zmm25, %zmm8
	vmovaps        160(%r12), %zmm5
	vfmadd213pd   %zmm2, %zmm24, %zmm14
	vmovnrngoaps  %zmm12,    0(%r13)
	vfmadd231pd   %zmm23, %zmm25, %zmm11
	addq          $-1, %r10
	vfmadd213pd   %zmm5, %zmm24, %zmm17
	vmovnrngoaps  %zmm13,   32(%r13)
	orq           %r10, %r10
	vmovnrngoaps  %zmm14,   64(%r13)
	vmovnrngoaps  %zmm15,   96(%r13)
	vmovnrngoaps  %zmm16,  128(%r13)
	vmovnrngoaps  %zmm17,  160(%r13)
	vmovnrngoaps  %zmm18,  192(%r13)
	vmovnrngoaps  %zmm19,  224(%r13)
	vmovnrngoaps  %zmm20,  256(%r13)
	vmovnrngoaps  %zmm21,  288(%r13)
	vmovnrngoaps  %zmm22,  320(%r13)
	vmovnrngoaps  %zmm23,  352(%r13)
	jg     _vmx_cayley_inv_s_lab1
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
	vmovnrngoaps  %zmm6,  192(%r13)
	vmulpd        %zmm11, %zmm26, %zmm11
	vmulpd        %zmm12, %zmm26, %zmm0
	vmovnrngoaps  %zmm7,  224(%r13)
	vmulpd        %zmm13, %zmm26, %zmm1
	vmulpd        %zmm14, %zmm26, %zmm2
	vmovnrngoaps  %zmm8,  256(%r13)
	vmulpd        %zmm15, %zmm26, %zmm3
	vmulpd        %zmm16, %zmm26, %zmm4
	vmovnrngoaps  %zmm0,    0(%r13)
	vmulpd        %zmm17, %zmm26, %zmm5
	vmovnrngoaps  %zmm1,   32(%r13)
	vmovnrngoaps  %zmm2,   64(%r13)
	vmovnrngoaps  %zmm3,   96(%r13)
	vmovnrngoaps  %zmm4,  128(%r13)
	vmovnrngoaps  %zmm5,  160(%r13)
	vmovnrngoaps  %zmm9,  288(%r13)
	vmovnrngoaps  %zmm10,  320(%r13)
	vmovnrngoaps  %zmm11,  352(%r13)
_vmx_cayley_inv_s_lab2:
	vbroadcastsd     0(%rbx), %zmm26
	vbroadcastsd     8(%rbx), %zmm24
	vbroadcastsd    16(%rbx), %zmm25
	addq          $24, %rbx
	addq          $-384, %r13
	addq          $-1, %r10
	vmovaps          0(%r13), %zmm12
	vmovaps        192(%r13), %zmm18
	vmovaps         96(%r13), %zmm15
	vmulpd        %zmm12, %zmm26, %zmm12
	vmovaps        288(%r13), %zmm21
	vmulpd        %zmm18, %zmm26, %zmm18
	vmovaps         32(%r13), %zmm13
	vmulpd        %zmm15, %zmm26, %zmm15
	vmovaps        224(%r13), %zmm19
	vmulpd        %zmm21, %zmm26, %zmm21
	vmovaps        128(%r13), %zmm16
	vmulpd        %zmm13, %zmm26, %zmm13
	vmovaps        320(%r13), %zmm22
	vmulpd        %zmm19, %zmm26, %zmm19
	vmovaps         64(%r13), %zmm14
	vmulpd        %zmm16, %zmm26, %zmm16
	vmovaps        256(%r13), %zmm20
	vmulpd        %zmm22, %zmm26, %zmm22
	vmovaps        160(%r13), %zmm17
	vmulpd        %zmm14, %zmm26, %zmm14
	vmovaps        352(%r13), %zmm23
	vmulpd        %zmm20, %zmm26, %zmm20
	orq           %r10, %r10
	vmulpd        %zmm17, %zmm26, %zmm17
	vmulpd        %zmm23, %zmm26, %zmm23
	vfmadd231pd   %zmm0, %zmm25, %zmm12
	vfmadd231pd   %zmm3, %zmm25, %zmm15
	vfmadd231pd   %zmm1, %zmm25, %zmm13
	vfmadd231pd   %zmm4, %zmm25, %zmm16
	vfmadd231pd   %zmm2, %zmm25, %zmm14
	vmovnrngoaps  %zmm12,    0(%r13)
	vfmadd231pd   %zmm5, %zmm25, %zmm17
	vfmadd213pd   %zmm18, %zmm24, %zmm6
	vmovnrngoaps  %zmm13,   32(%r13)
	vfmadd213pd   %zmm21, %zmm24, %zmm9
	vfmadd213pd   %zmm19, %zmm24, %zmm7
	vmovnrngoaps  %zmm14,   64(%r13)
	vfmadd213pd   %zmm22, %zmm24, %zmm10
	vfmadd213pd   %zmm20, %zmm24, %zmm8
	vmovnrngoaps  %zmm6,  192(%r13)
	vfmadd213pd   %zmm23, %zmm24, %zmm11
	vmovnrngoaps  %zmm7,  224(%r13)
	vmovnrngoaps  %zmm8,  256(%r13)
	vmovnrngoaps  %zmm9,  288(%r13)
	vmovnrngoaps  %zmm10,  320(%r13)
	vmovnrngoaps  %zmm11,  352(%r13)
	vmovnrngoaps  %zmm15,   96(%r13)
	vmovnrngoaps  %zmm16,  128(%r13)
	vmovnrngoaps  %zmm17,  160(%r13)
	jg     _vmx_cayley_inv_s_lab2
	imulq         $384, %rbp, %r12
	addq          $-1, %r9
	addq          %r12, %r8
	addq          %r12, %rcx
	orq           %r9, %r9
	jg     _vmx_cayley_inv_s_lab0
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
.type vmx_cayley_inv_s,@function
.size vmx_cayley_inv_s,.-vmx_cayley_inv_s
.data
