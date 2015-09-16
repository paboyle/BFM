/* 
 * BAGEL machine generated output.  
 *   
 * It is provided under the GNU pubic License V2  
 * It is provided as is, and is not guaranteed fit for any purpose.
 * BAGEL was written by Peter Boyle  
 */  
.text
.align 16,0x90
.globl vmx_vaxpby_s
vmx_vaxpby_s:
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
	movq          %rdi, %r11
	movq          %rsi, %r13
	movq          %rdx, %r14
	movq          %rcx, %rbp
	movq          %r8, %r10
	movq          %r9, %r12
	vbroadcastsd     0(%r13), %zmm0
	vbroadcastsd     0(%r14), %zmm1
	orq           %r12, %r12
	jle    _vmx_vaxpby_s_lab0
_vmx_vaxpby_s_lab1:
	movq          %r11, %r15
	vmovaps          0(%rbp), %zmm2
	vmovaps          0(%r10), %zmm8
	vmovaps         64(%rbp), %zmm4
	vmovaps         64(%r10), %zmm10
	vmovaps        128(%rbp), %zmm6
	vmovaps        128(%r10), %zmm12
	vmovaps         32(%rbp), %zmm3
	vmovaps         32(%r10), %zmm9
	vmovaps         96(%rbp), %zmm5
	vmovaps         96(%r10), %zmm11
	vmovaps        160(%rbp), %zmm7
	vmovaps        160(%r10), %zmm13
	vmulpd        %zmm8, %zmm1, %zmm8
	vmulpd        %zmm9, %zmm1, %zmm9
	vmulpd        %zmm10, %zmm1, %zmm10
	vmulpd        %zmm11, %zmm1, %zmm11
	vmulpd        %zmm12, %zmm1, %zmm12
	vmulpd        %zmm13, %zmm1, %zmm13
	vmovapd       %zmm8, %zmm14
	vfmadd231pd   %zmm2, %zmm0, %zmm14
	vmovapd       %zmm9, %zmm15
	vfmadd231pd   %zmm3, %zmm0, %zmm15
	vmovapd       %zmm10, %zmm16
	vfmadd231pd   %zmm4, %zmm0, %zmm16
	vmovapd       %zmm11, %zmm17
	vfmadd231pd   %zmm5, %zmm0, %zmm17
	vmovapd       %zmm12, %zmm18
	vfmadd231pd   %zmm6, %zmm0, %zmm18
	vmovapd       %zmm13, %zmm19
	vfmadd231pd   %zmm7, %zmm0, %zmm19
	addq          $192, %rbp
	addq          $192, %r10
	vmovnrngoaps  %zmm14,    0(%r15)
	vmovnrngoaps  %zmm15,   32(%r15)
	vmovnrngoaps  %zmm16,   64(%r15)
	vmovnrngoaps  %zmm17,   96(%r15)
	vmovnrngoaps  %zmm18,  128(%r15)
	vmovnrngoaps  %zmm19,  160(%r15)
	addq          $192, %r11
	addq          $-1, %r12
	orq           %r12, %r12
	jg     _vmx_vaxpby_s_lab1
_vmx_vaxpby_s_lab0:
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
.type vmx_vaxpby_s,@function
.size vmx_vaxpby_s,.-vmx_vaxpby_s
.data
