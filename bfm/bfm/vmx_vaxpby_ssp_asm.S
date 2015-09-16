/* 
 * BAGEL machine generated output.  
 *   
 * It is provided under the GNU pubic License V2  
 * It is provided as is, and is not guaranteed fit for any purpose.
 * BAGEL was written by Peter Boyle  
 */  
.text
.align 16,0x90
.globl vmx_vaxpby_ssp
vmx_vaxpby_ssp:
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
	movq          %rdi, %rcx
	movq             0(%rcx), %r10
	addq          $8, %rcx
	movq             0(%rcx), %r8
	addq          $8, %rcx
	movq             0(%rcx), %r9
	addq          $8, %rcx
	movq             0(%rcx), %rax
	addq          $8, %rcx
	movq             0(%rcx), %rbx
	addq          $8, %rcx
	movq             0(%rcx), %rbp
	vbroadcastsd     0(%rbp), %zmm0
	addq          $8, %rbp
	vbroadcastsd     0(%rbp), %zmm2
	addq          $8, %rbp
	vbroadcastsd     0(%rbp), %zmm1
	addq          $8, %rbp
	vbroadcastsd     0(%rbp), %zmm3
	addq          $8, %rbp
	orq           %rax, %rax
	jle    _vmx_vaxpby_ssp_lab0
_vmx_vaxpby_ssp_lab1:
	vmovapd          0(%r8), %zmm4
	vmovapd          0(%r9), %zmm10
	vmovapd        128(%r8), %zmm6
	vmovapd        128(%r9), %zmm12
	vmovapd        256(%r8), %zmm8
	vmovapd        256(%r9), %zmm14
	vmovapd         64(%r8), %zmm5
	vmovapd         64(%r9), %zmm11
	vmovapd        192(%r8), %zmm7
	vmovapd        192(%r9), %zmm13
	vmovapd        320(%r8), %zmm9
	vmovapd        320(%r9), %zmm15
	vmulpd        %zmm10, %zmm1, %zmm10
	vmovapd       %zmm10, %zmm16
	vfmadd231pd   %zmm4, %zmm0, %zmm16
	vmovnrngoapd  %zmm16,    0(%r10)
	vmulpd        %zmm11, %zmm1, %zmm11
	vmovapd       %zmm11, %zmm17
	vfmadd231pd   %zmm5, %zmm0, %zmm17
	vmovnrngoapd  %zmm17,   64(%r10)
	vmulpd        %zmm12, %zmm1, %zmm12
	vmovapd       %zmm12, %zmm18
	vfmadd231pd   %zmm6, %zmm0, %zmm18
	vmovnrngoapd  %zmm18,  128(%r10)
	vmulpd        %zmm13, %zmm1, %zmm13
	vmovapd       %zmm13, %zmm19
	vfmadd231pd   %zmm7, %zmm0, %zmm19
	vmovnrngoapd  %zmm19,  192(%r10)
	vmulpd        %zmm14, %zmm1, %zmm14
	vmovapd       %zmm14, %zmm20
	vfmadd231pd   %zmm8, %zmm0, %zmm20
	vmovnrngoapd  %zmm20,  256(%r10)
	vmulpd        %zmm15, %zmm1, %zmm15
	vmovapd       %zmm15, %zmm21
	vfmadd231pd   %zmm9, %zmm0, %zmm21
	vmovnrngoapd  %zmm21,  320(%r10)
	vmovapd        384(%r8), %zmm4
	vmovapd        384(%r9), %zmm10
	vmovapd        512(%r8), %zmm6
	vmovapd        512(%r9), %zmm12
	vmovapd        640(%r8), %zmm8
	vmovapd        640(%r9), %zmm14
	vmovapd        448(%r8), %zmm5
	vmovapd        448(%r9), %zmm11
	vmovapd        576(%r8), %zmm7
	vmovapd        576(%r9), %zmm13
	vmovapd        704(%r8), %zmm9
	vmovapd        704(%r9), %zmm15
	vmulpd        %zmm10, %zmm3, %zmm10
	vmovapd       %zmm10, %zmm16
	vfmadd231pd   %zmm4, %zmm2, %zmm16
	vmovnrngoapd  %zmm16,  384(%r10)
	vmulpd        %zmm11, %zmm3, %zmm11
	vmovapd       %zmm11, %zmm17
	vfmadd231pd   %zmm5, %zmm2, %zmm17
	vmovnrngoapd  %zmm17,  448(%r10)
	vmulpd        %zmm12, %zmm3, %zmm12
	vmovapd       %zmm12, %zmm18
	vfmadd231pd   %zmm6, %zmm2, %zmm18
	vmovnrngoapd  %zmm18,  512(%r10)
	vmulpd        %zmm13, %zmm3, %zmm13
	vmovapd       %zmm13, %zmm19
	vfmadd231pd   %zmm7, %zmm2, %zmm19
	vmovnrngoapd  %zmm19,  576(%r10)
	vmulpd        %zmm14, %zmm3, %zmm14
	vmovapd       %zmm14, %zmm20
	vfmadd231pd   %zmm8, %zmm2, %zmm20
	vmovnrngoapd  %zmm20,  640(%r10)
	vmulpd        %zmm15, %zmm3, %zmm15
	vmovapd       %zmm15, %zmm21
	vfmadd231pd   %zmm9, %zmm2, %zmm21
	vmovnrngoapd  %zmm21,  704(%r10)
	addq          %rbx, %r8
	addq          %rbx, %r9
	addq          %rbx, %r10
	addq          $-1, %rax
	orq           %rax, %rax
	jg     _vmx_vaxpby_ssp_lab1
_vmx_vaxpby_ssp_lab0:
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
.type vmx_vaxpby_ssp,@function
.size vmx_vaxpby_ssp,.-vmx_vaxpby_ssp
.data
