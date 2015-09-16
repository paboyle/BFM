/* 
 * BAGEL machine generated output.  
 *   
 * It is provided under the GNU pubic License V2  
 * It is provided as is, and is not guaranteed fit for any purpose.
 * BAGEL was written by Peter Boyle  
 */  
.text
.align 16,0x90
.globl vmx_caxpy_s
vmx_caxpy_s:
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
	movq          %rdx, %rbp
	movq          %rcx, %r10
	movq          %r8, %r12
	vbroadcastsd        0(%r13), %zmm0{%k7}
	vbroadcastsd        8(%r13), %zmm0{%k6}
	orq           %r12, %r12
	jle    _vmx_caxpy_s_lab0
_vmx_caxpy_s_lab1:
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
	vmulpd        %zmm2, %zmm0, %zmm20
	vmulpd        %zmm2{cdab}, %zmm0, %zmm21
	vsubpd        %zmm20{cdab}, %zmm20, %zmm20{%k7}
	vaddpd        %zmm21{cdab}, %zmm21, %zmm20{%k6}
	vaddpd        %zmm8, %zmm20, %zmm14
	vmulpd        %zmm3, %zmm0, %zmm22
	vmulpd        %zmm3{cdab}, %zmm0, %zmm23
	vsubpd        %zmm22{cdab}, %zmm22, %zmm22{%k7}
	vaddpd        %zmm23{cdab}, %zmm23, %zmm22{%k6}
	vaddpd        %zmm9, %zmm22, %zmm15
	vmulpd        %zmm4, %zmm0, %zmm24
	vmulpd        %zmm4{cdab}, %zmm0, %zmm25
	vsubpd        %zmm24{cdab}, %zmm24, %zmm24{%k7}
	vaddpd        %zmm25{cdab}, %zmm25, %zmm24{%k6}
	vaddpd        %zmm10, %zmm24, %zmm16
	vmulpd        %zmm5, %zmm0, %zmm26
	vmulpd        %zmm5{cdab}, %zmm0, %zmm27
	vsubpd        %zmm26{cdab}, %zmm26, %zmm26{%k7}
	vaddpd        %zmm27{cdab}, %zmm27, %zmm26{%k6}
	vaddpd        %zmm11, %zmm26, %zmm17
	vmulpd        %zmm6, %zmm0, %zmm28
	vmulpd        %zmm6{cdab}, %zmm0, %zmm29
	vsubpd        %zmm28{cdab}, %zmm28, %zmm28{%k7}
	vaddpd        %zmm29{cdab}, %zmm29, %zmm28{%k6}
	vaddpd        %zmm12, %zmm28, %zmm18
	vmulpd        %zmm7, %zmm0, %zmm30
	vmulpd        %zmm7{cdab}, %zmm0, %zmm31
	vsubpd        %zmm30{cdab}, %zmm30, %zmm30{%k7}
	vaddpd        %zmm31{cdab}, %zmm31, %zmm30{%k6}
	vaddpd        %zmm13, %zmm30, %zmm19
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
	jg     _vmx_caxpy_s_lab1
_vmx_caxpy_s_lab0:
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
.type vmx_caxpy_s,@function
.size vmx_caxpy_s,.-vmx_caxpy_s
.data
