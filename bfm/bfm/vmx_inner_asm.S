/* 
 * BAGEL machine generated output.  
 *   
 * It is provided under the GNU pubic License V2  
 * It is provided as is, and is not guaranteed fit for any purpose.
 * BAGEL was written by Peter Boyle  
 */  
.text
.align 16,0x90
.globl vmx_inner
vmx_inner:
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
	movq          %rdi, %rbx
	movq          %rsi, %rbp
	movq          %rdx, %r10
	movq          %rcx, %r11
_vmx_inner_lab0:
	vmovapd          0(%r11), %zmm12
	vmovapd          0(%r11), %zmm13
	vmovapd          0(%r11), %zmm14
	vmovapd          0(%r11), %zmm15
	vmovapd          0(%r11), %zmm16
	vmovapd          0(%r11), %zmm17
	vmovapd          0(%r11), %zmm18
	vbroadcastsd     0(%r11), %zmm20
	vbroadcastsd     0(%r11), %zmm21
	orq           %r10, %r10
	jle    _vmx_inner_lab1
_vmx_inner_lab2:
	vmovapd          0(%rbx), %zmm0
	vmovapd          0(%rbp), %zmm6
	vmovapd        128(%rbx), %zmm2
	vmovapd        128(%rbp), %zmm8
	vmovapd        256(%rbx), %zmm4
	vmovapd        256(%rbp), %zmm10
	vmovapd         64(%rbx), %zmm1
	vmovapd         64(%rbp), %zmm7
	vmovapd        192(%rbx), %zmm3
	vmovapd        192(%rbp), %zmm9
	vmovapd        320(%rbx), %zmm5
	vmovapd        320(%rbp), %zmm11
	addq          $384, %rbx
	addq          $384, %rbp
	vmulpd        %zmm6, %zmm0, %zmm24
	vmulpd        %zmm6{cdab}, %zmm0, %zmm25
	vaddpd        %zmm24{cdab}, %zmm24, %zmm24{%k7}
	vsubrpd       %zmm25{cdab}, %zmm25, %zmm24{%k6}
	vaddpd        %zmm12, %zmm24, %zmm12
	vmulpd        %zmm7, %zmm1, %zmm26
	vmulpd        %zmm7{cdab}, %zmm1, %zmm27
	vaddpd        %zmm26{cdab}, %zmm26, %zmm26{%k7}
	vsubrpd       %zmm27{cdab}, %zmm27, %zmm26{%k6}
	vaddpd        %zmm13, %zmm26, %zmm13
	vmulpd        %zmm8, %zmm2, %zmm28
	vmulpd        %zmm8{cdab}, %zmm2, %zmm29
	vaddpd        %zmm28{cdab}, %zmm28, %zmm28{%k7}
	vsubrpd       %zmm29{cdab}, %zmm29, %zmm28{%k6}
	vaddpd        %zmm14, %zmm28, %zmm14
	vmulpd        %zmm9, %zmm3, %zmm30
	vmulpd        %zmm9{cdab}, %zmm3, %zmm31
	vaddpd        %zmm30{cdab}, %zmm30, %zmm30{%k7}
	vsubrpd       %zmm31{cdab}, %zmm31, %zmm30{%k6}
	vaddpd        %zmm15, %zmm30, %zmm15
	vmulpd        %zmm10, %zmm4, %zmm24
	vmulpd        %zmm10{cdab}, %zmm4, %zmm25
	vaddpd        %zmm24{cdab}, %zmm24, %zmm24{%k7}
	vsubrpd       %zmm25{cdab}, %zmm25, %zmm24{%k6}
	vaddpd        %zmm16, %zmm24, %zmm16
	vmulpd        %zmm11, %zmm5, %zmm24
	vmulpd        %zmm11{cdab}, %zmm5, %zmm25
	vaddpd        %zmm24{cdab}, %zmm24, %zmm24{%k7}
	vsubrpd       %zmm25{cdab}, %zmm25, %zmm24{%k6}
	vaddpd        %zmm17, %zmm24, %zmm17
	addq          $-1, %r10
	orq           %r10, %r10
	jg     _vmx_inner_lab2
	vaddpd        %zmm13, %zmm12, %zmm12
	vaddpd        %zmm14, %zmm12, %zmm12
	vaddpd        %zmm15, %zmm12, %zmm12
	vaddpd        %zmm16, %zmm12, %zmm12
	vaddpd        %zmm17, %zmm12, %zmm12
	vmovnrngoapd  %zmm12,    0(%r11)
	vbroadcastsd     0(%r11), %zmm22
	vaddpd        %zmm22, %zmm20, %zmm20{%k5}
	vbroadcastsd     8(%r11), %zmm23
	vaddpd        %zmm23, %zmm21, %zmm21{%k5}
	vbroadcastsd    16(%r11), %zmm22
	vaddpd        %zmm22, %zmm20, %zmm20{%k5}
	vbroadcastsd    24(%r11), %zmm23
	vaddpd        %zmm23, %zmm21, %zmm21{%k5}
	vbroadcastsd    32(%r11), %zmm22
	vaddpd        %zmm22, %zmm20, %zmm20{%k5}
	vbroadcastsd    40(%r11), %zmm23
	vaddpd        %zmm23, %zmm21, %zmm21{%k5}
	vbroadcastsd    48(%r11), %zmm22
	vaddpd        %zmm22, %zmm20, %zmm20{%k5}
	vbroadcastsd    56(%r11), %zmm23
	vaddpd        %zmm23, %zmm21, %zmm21{%k5}
	vpackstorelpd %zmm20,    0(%r11){%k5}
	vpackstorelpd %zmm21,    8(%r11){%k5}
_vmx_inner_lab1:
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
.type vmx_inner,@function
.size vmx_inner,.-vmx_inner
.data
