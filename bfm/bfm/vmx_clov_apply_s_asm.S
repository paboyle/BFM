/* 
 * BAGEL machine generated output.  
 *   
 * It is provided under the GNU pubic License V2  
 * It is provided as is, and is not guaranteed fit for any purpose.
 * BAGEL was written by Peter Boyle  
 */  
.text
.align 16,0x90
.globl vmx_clov_apply_s
vmx_clov_apply_s:
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
	movq             0(%rbx), %r8
	addq          $8, %rbx
	movq             0(%rbx), %r9
	addq          $8, %rbx
	movq             0(%rbx), %rcx
	addq          $8, %rbx
	movq             0(%rbx), %rax
	orq           %rax, %rax
	jle    _vmx_clov_apply_s_lab0
_vmx_clov_apply_s_lab1:
	vmovaps          0(%r9), %zmm6
	vmovaps         32(%r9), %zmm7
	vmovaps         64(%r9), %zmm8
	vmovaps         96(%r9), %zmm9
	vmovaps        128(%r9), %zmm10
	vmovaps        160(%r9), %zmm11
	vmovaps          0(%rcx), %zmm12
	vmovaps         32(%rcx), %zmm13
	vmovaps         64(%rcx), %zmm14
	vmulpd        %zmm6, %zmm12, %zmm0
	vmovaps         96(%rcx), %zmm15
	vmulpd        %zmm6{cdab}, %zmm12, %zmm27
	vmovaps        128(%rcx), %zmm16
	vmulpd        %zmm7, %zmm13, %zmm1
	vmovaps        160(%rcx), %zmm17
	vmulpd        %zmm7{cdab}, %zmm13, %zmm28
	vmovaps        192(%rcx), %zmm12
	vsubpd        %zmm0{cdab}, %zmm0, %zmm0{%k7}
	vmovaps        224(%rcx), %zmm13
	vmulpd        %zmm8, %zmm14, %zmm2
	vmovaps        384(%rcx), %zmm18
	vsubpd        %zmm1{cdab}, %zmm1, %zmm1{%k7}
	vmovaps        416(%rcx), %zmm19
	vaddpd        %zmm27{cdab}, %zmm27, %zmm0{%k6}
	vmovaps        448(%rcx), %zmm20
	vmulpd        %zmm8{cdab}, %zmm14, %zmm29
	vmovaps        480(%rcx), %zmm21
	vaddpd        %zmm28{cdab}, %zmm28, %zmm1{%k6}
	vmovaps        256(%rcx), %zmm14
	vsubpd        %zmm2{cdab}, %zmm2, %zmm2{%k7}
	vmovaps        512(%rcx), %zmm22
	vmulpd        %zmm9, %zmm15, %zmm3
	vmovaps        544(%rcx), %zmm23
	vmulpd        %zmm9{cdab}, %zmm15, %zmm30
	vmovaps        576(%rcx), %zmm24
	vmulpd        %zmm10, %zmm16, %zmm4
	vmovaps        288(%rcx), %zmm15
	vaddpd        %zmm29{cdab}, %zmm29, %zmm2{%k6}
	vmovaps        608(%rcx), %zmm25
	vsubpd        %zmm3{cdab}, %zmm3, %zmm3{%k7}
	vmovaps        640(%rcx), %zmm26
	vmulpd        %zmm10{cdab}, %zmm16, %zmm31
	addq          $192, %r9
	vsubpd        %zmm4{cdab}, %zmm4, %zmm4{%k7}
	vmovaps        320(%rcx), %zmm16
	vaddpd        %zmm30{cdab}, %zmm30, %zmm3{%k6}
	addq          $-1, %rax
	vmulpd        %zmm11, %zmm17, %zmm5
	vmulpd        %zmm11{cdab}, %zmm17, %zmm27
	vmulpd        %zmm7{cdab}, %zmm12, %zmm28
	vmovaps        352(%rcx), %zmm17
	vaddpd        %zmm31{cdab}, %zmm31, %zmm4{%k6}
	addq          $672, %rcx
	vsubpd        %zmm5{cdab}, %zmm5, %zmm5{%k7}
	vaddpd        %zmm27{cdab}, %zmm27, %zmm5{%k6}
	vmulpd        %zmm7, %zmm12, %zmm27
	vaddpd        %zmm27{cdab}, %zmm27, %zmm27{%k7}
	vsubrpd       %zmm28{cdab}, %zmm28, %zmm27{%k6}
	vmulpd        %zmm8{cdab}, %zmm13, %zmm28
	vaddpd        %zmm0, %zmm27, %zmm0
	vmulpd        %zmm8, %zmm13, %zmm27
	vaddpd        %zmm27{cdab}, %zmm27, %zmm27{%k7}
	vsubrpd       %zmm28{cdab}, %zmm28, %zmm27{%k6}
	vmulpd        %zmm9{cdab}, %zmm15, %zmm28
	vaddpd        %zmm0, %zmm27, %zmm0
	vmulpd        %zmm9, %zmm15, %zmm27
	vaddpd        %zmm27{cdab}, %zmm27, %zmm27{%k7}
	vsubrpd       %zmm28{cdab}, %zmm28, %zmm27{%k6}
	vmulpd        %zmm10{cdab}, %zmm18, %zmm28
	vaddpd        %zmm0, %zmm27, %zmm0
	vmulpd        %zmm10, %zmm18, %zmm27
	vaddpd        %zmm27{cdab}, %zmm27, %zmm27{%k7}
	vsubrpd       %zmm28{cdab}, %zmm28, %zmm27{%k6}
	vmulpd        %zmm11{cdab}, %zmm22, %zmm28
	vaddpd        %zmm0, %zmm27, %zmm0
	vmulpd        %zmm11, %zmm22, %zmm27
	vaddpd        %zmm27{cdab}, %zmm27, %zmm27{%k7}
	vsubrpd       %zmm28{cdab}, %zmm28, %zmm27{%k6}
	vmulpd        %zmm6{cdab}, %zmm12, %zmm28
	vaddpd        %zmm0, %zmm27, %zmm0
	vmulpd        %zmm6, %zmm12, %zmm27
	vmovnrngoaps  %zmm0,    0(%r8)
	vsubpd        %zmm27{cdab}, %zmm27, %zmm27{%k7}
	vaddpd        %zmm28{cdab}, %zmm28, %zmm27{%k6}
	vmulpd        %zmm8{cdab}, %zmm14, %zmm28
	vaddpd        %zmm1, %zmm27, %zmm1
	vmulpd        %zmm8, %zmm14, %zmm27
	vaddpd        %zmm27{cdab}, %zmm27, %zmm27{%k7}
	vsubrpd       %zmm28{cdab}, %zmm28, %zmm27{%k6}
	vmulpd        %zmm9{cdab}, %zmm16, %zmm28
	vaddpd        %zmm1, %zmm27, %zmm1
	vmulpd        %zmm9, %zmm16, %zmm27
	vaddpd        %zmm27{cdab}, %zmm27, %zmm27{%k7}
	vsubrpd       %zmm28{cdab}, %zmm28, %zmm27{%k6}
	vmulpd        %zmm10{cdab}, %zmm19, %zmm28
	vaddpd        %zmm1, %zmm27, %zmm1
	vmulpd        %zmm10, %zmm19, %zmm27
	vaddpd        %zmm27{cdab}, %zmm27, %zmm27{%k7}
	vsubrpd       %zmm28{cdab}, %zmm28, %zmm27{%k6}
	vmulpd        %zmm11{cdab}, %zmm23, %zmm28
	vaddpd        %zmm1, %zmm27, %zmm1
	vmulpd        %zmm11, %zmm23, %zmm27
	vaddpd        %zmm27{cdab}, %zmm27, %zmm27{%k7}
	vsubrpd       %zmm28{cdab}, %zmm28, %zmm27{%k6}
	vmulpd        %zmm6{cdab}, %zmm13, %zmm28
	vaddpd        %zmm1, %zmm27, %zmm1
	vmulpd        %zmm6, %zmm13, %zmm27
	vmovnrngoaps  %zmm1,   32(%r8)
	vsubpd        %zmm27{cdab}, %zmm27, %zmm27{%k7}
	vaddpd        %zmm28{cdab}, %zmm28, %zmm27{%k6}
	vmulpd        %zmm7{cdab}, %zmm14, %zmm28
	vaddpd        %zmm2, %zmm27, %zmm2
	vmulpd        %zmm7, %zmm14, %zmm27
	vsubpd        %zmm27{cdab}, %zmm27, %zmm27{%k7}
	vaddpd        %zmm28{cdab}, %zmm28, %zmm27{%k6}
	vmulpd        %zmm9{cdab}, %zmm17, %zmm28
	vaddpd        %zmm2, %zmm27, %zmm2
	vmulpd        %zmm9, %zmm17, %zmm27
	vaddpd        %zmm27{cdab}, %zmm27, %zmm27{%k7}
	vsubrpd       %zmm28{cdab}, %zmm28, %zmm27{%k6}
	vmulpd        %zmm10{cdab}, %zmm20, %zmm28
	vaddpd        %zmm2, %zmm27, %zmm2
	vmulpd        %zmm10, %zmm20, %zmm27
	vaddpd        %zmm27{cdab}, %zmm27, %zmm27{%k7}
	vsubrpd       %zmm28{cdab}, %zmm28, %zmm27{%k6}
	vmulpd        %zmm11{cdab}, %zmm24, %zmm28
	vaddpd        %zmm2, %zmm27, %zmm2
	vmulpd        %zmm11, %zmm24, %zmm27
	vaddpd        %zmm27{cdab}, %zmm27, %zmm27{%k7}
	vsubrpd       %zmm28{cdab}, %zmm28, %zmm27{%k6}
	vmulpd        %zmm6{cdab}, %zmm15, %zmm28
	vaddpd        %zmm2, %zmm27, %zmm2
	vmulpd        %zmm6, %zmm15, %zmm27
	vmovnrngoaps  %zmm2,   64(%r8)
	vsubpd        %zmm27{cdab}, %zmm27, %zmm27{%k7}
	vaddpd        %zmm28{cdab}, %zmm28, %zmm27{%k6}
	vmulpd        %zmm7{cdab}, %zmm16, %zmm28
	vaddpd        %zmm3, %zmm27, %zmm3
	vmulpd        %zmm7, %zmm16, %zmm27
	vsubpd        %zmm27{cdab}, %zmm27, %zmm27{%k7}
	vaddpd        %zmm28{cdab}, %zmm28, %zmm27{%k6}
	vmulpd        %zmm8{cdab}, %zmm17, %zmm28
	vaddpd        %zmm3, %zmm27, %zmm3
	vmulpd        %zmm8, %zmm17, %zmm27
	vsubpd        %zmm27{cdab}, %zmm27, %zmm27{%k7}
	vaddpd        %zmm28{cdab}, %zmm28, %zmm27{%k6}
	vmulpd        %zmm10{cdab}, %zmm21, %zmm28
	vaddpd        %zmm3, %zmm27, %zmm3
	vmulpd        %zmm10, %zmm21, %zmm27
	vaddpd        %zmm27{cdab}, %zmm27, %zmm27{%k7}
	vsubrpd       %zmm28{cdab}, %zmm28, %zmm27{%k6}
	vmulpd        %zmm11{cdab}, %zmm25, %zmm28
	vaddpd        %zmm3, %zmm27, %zmm3
	vmulpd        %zmm11, %zmm25, %zmm27
	vaddpd        %zmm27{cdab}, %zmm27, %zmm27{%k7}
	vsubrpd       %zmm28{cdab}, %zmm28, %zmm27{%k6}
	vmulpd        %zmm6{cdab}, %zmm18, %zmm28
	vaddpd        %zmm3, %zmm27, %zmm3
	vmulpd        %zmm6, %zmm18, %zmm27
	vmovnrngoaps  %zmm3,   96(%r8)
	vsubpd        %zmm27{cdab}, %zmm27, %zmm27{%k7}
	vaddpd        %zmm28{cdab}, %zmm28, %zmm27{%k6}
	vmulpd        %zmm7{cdab}, %zmm19, %zmm28
	vaddpd        %zmm4, %zmm27, %zmm4
	vmulpd        %zmm7, %zmm19, %zmm27
	vsubpd        %zmm27{cdab}, %zmm27, %zmm27{%k7}
	vaddpd        %zmm28{cdab}, %zmm28, %zmm27{%k6}
	vmulpd        %zmm8{cdab}, %zmm20, %zmm28
	vaddpd        %zmm4, %zmm27, %zmm4
	vmulpd        %zmm8, %zmm20, %zmm27
	vsubpd        %zmm27{cdab}, %zmm27, %zmm27{%k7}
	vaddpd        %zmm28{cdab}, %zmm28, %zmm27{%k6}
	vmulpd        %zmm9{cdab}, %zmm21, %zmm28
	vaddpd        %zmm4, %zmm27, %zmm4
	vmulpd        %zmm9, %zmm21, %zmm27
	vsubpd        %zmm27{cdab}, %zmm27, %zmm27{%k7}
	vaddpd        %zmm28{cdab}, %zmm28, %zmm27{%k6}
	vmulpd        %zmm11{cdab}, %zmm26, %zmm28
	vaddpd        %zmm4, %zmm27, %zmm4
	vmulpd        %zmm11, %zmm26, %zmm27
	vaddpd        %zmm27{cdab}, %zmm27, %zmm27{%k7}
	vsubrpd       %zmm28{cdab}, %zmm28, %zmm27{%k6}
	vmulpd        %zmm6{cdab}, %zmm22, %zmm28
	vaddpd        %zmm4, %zmm27, %zmm4
	vmulpd        %zmm6, %zmm22, %zmm27
	vmovnrngoaps  %zmm4,  128(%r8)
	vsubpd        %zmm27{cdab}, %zmm27, %zmm27{%k7}
	vaddpd        %zmm28{cdab}, %zmm28, %zmm27{%k6}
	vmulpd        %zmm7{cdab}, %zmm23, %zmm28
	vaddpd        %zmm5, %zmm27, %zmm5
	vmulpd        %zmm7, %zmm23, %zmm27
	vsubpd        %zmm27{cdab}, %zmm27, %zmm27{%k7}
	vaddpd        %zmm28{cdab}, %zmm28, %zmm27{%k6}
	vmulpd        %zmm8{cdab}, %zmm24, %zmm28
	vaddpd        %zmm5, %zmm27, %zmm5
	vmulpd        %zmm8, %zmm24, %zmm27
	vsubpd        %zmm27{cdab}, %zmm27, %zmm27{%k7}
	vaddpd        %zmm28{cdab}, %zmm28, %zmm27{%k6}
	vmulpd        %zmm9{cdab}, %zmm25, %zmm28
	vaddpd        %zmm5, %zmm27, %zmm5
	vmulpd        %zmm9, %zmm25, %zmm27
	vsubpd        %zmm27{cdab}, %zmm27, %zmm27{%k7}
	vaddpd        %zmm28{cdab}, %zmm28, %zmm27{%k6}
	vmulpd        %zmm10{cdab}, %zmm26, %zmm28
	vaddpd        %zmm5, %zmm27, %zmm5
	vmulpd        %zmm10, %zmm26, %zmm27
	vsubpd        %zmm27{cdab}, %zmm27, %zmm27{%k7}
	vaddpd        %zmm28{cdab}, %zmm28, %zmm27{%k6}
	vaddpd        %zmm5, %zmm27, %zmm5
	vmovnrngoaps  %zmm5,  160(%r8)
	addq          $192, %r8
	orq           %rax, %rax
	jg     _vmx_clov_apply_s_lab1
_vmx_clov_apply_s_lab0:
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
.type vmx_clov_apply_s,@function
.size vmx_clov_apply_s,.-vmx_clov_apply_s
.data
