/* 
 * BAGEL machine generated output.  
 *   
 * It is provided under the GNU pubic License V2  
 * It is provided as is, and is not guaranteed fit for any purpose.
 * BAGEL was written by Peter Boyle  
 */  
.text
.align 16,0x90
.globl vmx_cayley_s
vmx_cayley_s:
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
	movq          %rbx, 2112(%rsp)
	movq          %rbp, 2120(%rsp)
	movq          %r12, 2144(%rsp)
	movq          %r13, 2152(%rsp)
	movq          %r14, 2160(%rsp)
	movq          %r15, 2168(%rsp)
	movq          %rdi, %r13
	movq             0(%r13), %r8
	addq          $8, %r13
	movq             0(%r13), %rcx
	addq          $8, %r13
	movq             0(%r13), %r9
	addq          $8, %r13
	movq             0(%r13), %r11
	addq          $8, %r13
	movq             0(%r13), %rax
	addq          $8, %r13
	movq             0(%r13), %rbp
	addq          $8, %r13
_vmx_cayley_s_lab0:
	movq          %r8, %r14
	movq          %rbp, %r10
	movq          %rax, %rbx
	movq          %r11, %r12
_vmx_cayley_s_lab1:
	movq             0(%r10), %r15
	movq             8(%r10), %rdi
	imulq         $384, %r15, %r15
	imulq         $384, %rdi, %rdi
	addq          %r8, %r15
	addq          %r8, %rdi
	vbroadcastsd     0(%rbx), %zmm24
	vbroadcastsd     8(%rbx), %zmm26
	vbroadcastsd    16(%rbx), %zmm25
	vmovaps          0(%rdi), %zmm0
	vmovaps        192(%r15), %zmm6
	vmovaps          0(%r14), %zmm12
	vmovaps        192(%r14), %zmm18
	vmovaps         96(%rdi), %zmm3
	vmulpd        %zmm12, %zmm26, %zmm12
	vmovaps        288(%r15), %zmm9
	vmulpd        %zmm18, %zmm26, %zmm18
	vmovaps         96(%r14), %zmm15
	vfmadd231pd   %zmm0, %zmm25, %zmm12
	vmovaps        288(%r14), %zmm21
	vfmadd231pd   %zmm6, %zmm24, %zmm18
	vmovaps         32(%rdi), %zmm1
	vmulpd        %zmm15, %zmm26, %zmm15
	vmovaps        224(%r15), %zmm7
	vmulpd        %zmm21, %zmm26, %zmm21
	vmovaps         32(%r14), %zmm13
	vfmadd231pd   %zmm3, %zmm25, %zmm15
	vmovaps        224(%r14), %zmm19
	vfmadd231pd   %zmm9, %zmm24, %zmm21
	vmovaps        128(%rdi), %zmm4
	vmulpd        %zmm13, %zmm26, %zmm13
	vmovaps        320(%r15), %zmm10
	vmulpd        %zmm19, %zmm26, %zmm19
	vmovaps        128(%r14), %zmm16
	vfmadd231pd   %zmm1, %zmm25, %zmm13
	vmovaps        320(%r14), %zmm22
	vfmadd231pd   %zmm7, %zmm24, %zmm19
	vmovaps         64(%rdi), %zmm2
	vmulpd        %zmm16, %zmm26, %zmm16
	vmovaps        256(%r15), %zmm8
	vmulpd        %zmm22, %zmm26, %zmm22
	vmovaps         64(%r14), %zmm14
	vfmadd231pd   %zmm4, %zmm25, %zmm16
	vmovaps        256(%r14), %zmm20
	vfmadd231pd   %zmm10, %zmm24, %zmm22
	vmovaps        160(%rdi), %zmm5
	vmulpd        %zmm14, %zmm26, %zmm14
	vmovaps        352(%r15), %zmm11
	vmulpd        %zmm20, %zmm26, %zmm20
	vmovaps        160(%r14), %zmm17
	vfmadd231pd   %zmm2, %zmm25, %zmm14
	vmovaps        352(%r14), %zmm23
	vfmadd231pd   %zmm8, %zmm24, %zmm20
	vmovnrngoaps  %zmm12,    0(%rcx)
	vmulpd        %zmm17, %zmm26, %zmm17
	addq          $384, %r14
	vmulpd        %zmm23, %zmm26, %zmm23
	vmovnrngoaps  %zmm13,   32(%rcx)
	vfmadd231pd   %zmm5, %zmm25, %zmm17
	addq          $16, %r10
	vfmadd231pd   %zmm11, %zmm24, %zmm23
	vmovnrngoaps  %zmm14,   64(%rcx)
	addq          $24, %rbx
	vmovnrngoaps  %zmm15,   96(%rcx)
	addq          $-1, %r12
	vmovnrngoaps  %zmm16,  128(%rcx)
	vmovnrngoaps  %zmm17,  160(%rcx)
	vmovnrngoaps  %zmm18,  192(%rcx)
	vmovnrngoaps  %zmm19,  224(%rcx)
	vmovnrngoaps  %zmm20,  256(%rcx)
	vmovnrngoaps  %zmm21,  288(%rcx)
	vmovnrngoaps  %zmm22,  320(%rcx)
	vmovnrngoaps  %zmm23,  352(%rcx)
	addq          $384, %rcx
	orq           %r12, %r12
	jg     _vmx_cayley_s_lab1
	imulq         $384, %r11, %r14
	addq          $-1, %r9
	addq          %r14, %r8
	orq           %r9, %r9
	jg     _vmx_cayley_s_lab0
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
.type vmx_cayley_s,@function
.size vmx_cayley_s,.-vmx_cayley_s
.data
