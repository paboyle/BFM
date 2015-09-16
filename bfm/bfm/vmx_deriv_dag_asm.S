/* 
 * BAGEL machine generated output.  
 *   
 * It is provided under the GNU pubic License V2  
 * It is provided as is, and is not guaranteed fit for any purpose.
 * BAGEL was written by Peter Boyle  
 */  
.text
.align 16,0x90
.globl vmx_deriv_dag
vmx_deriv_dag:
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
	movq             0(%r13), %r9
	addq          $8, %r13
	movq             0(%r13), %r8
	addq          $8, %r13
	movq             0(%r13), %rsi
	addq          $8, %r13
	movq             0(%r13), %rax
	addq          $8, %r13
	movq             0(%r13), %r10
	addq          $8, %r13
	movq             0(%r13), %rbx
	addq          $8, %r13
	movq             0(%r13), %rbp
	addq          $8, %r13
	addq          $8, %r13
	vmovapd          0(%rbp), %zmm28
	movq             0(%r13), %r12
	addq          $8, %r13
	vmovapd        128(%rbp), %zmm27
	orq           %rax, %rax
	jle    _vmx_deriv_dag_lab0
_vmx_deriv_dag_lab1:
	movq            80(%rbx), %r15
	movq          %r10, %r11
	movq          $0,%r14
_vmx_deriv_dag_lab2:
	movq          %rsi, %rdx
	movq          %r15, %rcx
	andq          $2, %rcx
	orq           %rcx, %rcx
	jle    _vmx_deriv_dag_lab3
	movq             8(%rbx), %rcx
	addq          %r12, %rcx
	addq          %r14, %rcx
	vmovapd          0(%rcx), %zmm0
	vmovapd         64(%rcx), %zmm1
	vmovapd        128(%rcx), %zmm2
	vmovapd        192(%rcx), %zmm3
	vmovapd        256(%rcx), %zmm4
	vmovapd        320(%rcx), %zmm5
	jmp    _vmx_deriv_dag_lab4
_vmx_deriv_dag_lab3:
	movq             8(%rbx), %rcx
	addq          %r8, %rcx
	addq          %r14, %rcx
	addq          %r14, %rcx
	vmovapd          0(%rcx), %zmm12
	vmovapd         64(%rcx), %zmm13
	vmovapd        128(%rcx), %zmm14
	vmovapd        192(%rcx), %zmm15
	vmovapd        256(%rcx), %zmm16
	vmovapd        320(%rcx), %zmm17
	vmovapd        384(%rcx), %zmm18
	vmovapd        448(%rcx), %zmm19
	vmovapd        512(%rcx), %zmm20
	vsubpd        %zmm18{cdab}, %zmm15, %zmm3{%k7}
	vmovapd        576(%rcx), %zmm21
	vaddpd        %zmm18{cdab}, %zmm15, %zmm3{%k6}
	vmovapd        640(%rcx), %zmm22
	vsubpd        %zmm19{cdab}, %zmm16, %zmm4{%k7}
	vmovapd        704(%rcx), %zmm23
	vsubpd        %zmm21{cdab}, %zmm12, %zmm0{%k7}
	vaddpd        %zmm21{cdab}, %zmm12, %zmm0{%k6}
	vsubpd        %zmm22{cdab}, %zmm13, %zmm1{%k7}
	vaddpd        %zmm22{cdab}, %zmm13, %zmm1{%k6}
	vsubpd        %zmm23{cdab}, %zmm14, %zmm2{%k7}
	vaddpd        %zmm23{cdab}, %zmm14, %zmm2{%k6}
	vaddpd        %zmm19{cdab}, %zmm16, %zmm4{%k6}
	vsubpd        %zmm20{cdab}, %zmm17, %zmm5{%k7}
	vaddpd        %zmm20{cdab}, %zmm17, %zmm5{%k6}
_vmx_deriv_dag_lab4:
	vaddpd        %zmm3{cdab}, %zmm27, %zmm6{%k7}
	vmovapd          0(%r9), %zmm12
	vsubpd        %zmm3{cdab}, %zmm27, %zmm6{%k6}
	vmovapd        192(%r9), %zmm15
	vaddpd        %zmm0{cdab}, %zmm27, %zmm9{%k7}
	vmovapd        384(%r9), %zmm18
	vsubpd        %zmm0{cdab}, %zmm27, %zmm9{%k6}
	vmovapd        576(%r9), %zmm21
	vaddpd        %zmm4{cdab}, %zmm27, %zmm7{%k7}
	vmovapd         64(%r9), %zmm13
	vsubpd        %zmm4{cdab}, %zmm27, %zmm7{%k6}
	vmovapd        256(%r9), %zmm16
	vaddpd        %zmm1{cdab}, %zmm27, %zmm10{%k7}
	vmovapd        448(%r9), %zmm19
	vsubpd        %zmm1{cdab}, %zmm27, %zmm10{%k6}
	vmovapd        640(%r9), %zmm22
	vaddpd        %zmm5{cdab}, %zmm27, %zmm8{%k7}
	vmovapd        128(%r9), %zmm14
	vsubpd        %zmm5{cdab}, %zmm27, %zmm8{%k6}
	vmovapd        320(%r9), %zmm17
	vaddpd        %zmm2{cdab}, %zmm27, %zmm11{%k7}
	vmovapd        512(%r9), %zmm20
	vsubpd        %zmm2{cdab}, %zmm27, %zmm11{%k6}
	vmovapd        704(%r9), %zmm23
	vmulpd        %zmm0, %zmm12, %zmm29
	vmovapd          0(%rdx), %zmm24
	vmulpd        %zmm0{cdab}, %zmm12, %zmm30
	vmovapd        192(%rdx), %zmm25
	vmulpd        %zmm1, %zmm12, %zmm31
	vmovapd        384(%rdx), %zmm26
	movq          %r15, %rcx
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	andq          $8, %rcx
	vaddpd        %zmm31{cdab}, %zmm31, %zmm31{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm2{cdab}, %zmm12, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm1{cdab}, %zmm12, %zmm29
	vsubrpd       %zmm29{cdab}, %zmm29, %zmm31{%k6}
	vmulpd        %zmm2, %zmm12, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vaddpd        %zmm25, %zmm31, %zmm25
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm3{cdab}, %zmm15, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm3, %zmm15, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm4{cdab}, %zmm15, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm4, %zmm15, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm5{cdab}, %zmm15, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm5, %zmm15, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm6{cdab}, %zmm18, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm6, %zmm18, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm7{cdab}, %zmm18, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm7, %zmm18, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm8{cdab}, %zmm18, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm8, %zmm18, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm9{cdab}, %zmm21, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm9, %zmm21, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm10{cdab}, %zmm21, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm10, %zmm21, %zmm29
	vmovnrngoapd  %zmm24,    0(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vmovapd         64(%rdx), %zmm24
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm11{cdab}, %zmm21, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm11, %zmm21, %zmm29
	vmovnrngoapd  %zmm25,  192(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vmovapd        256(%rdx), %zmm25
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm0{cdab}, %zmm13, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm0, %zmm13, %zmm29
	vmovnrngoapd  %zmm26,  384(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vmovapd        448(%rdx), %zmm26
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm1{cdab}, %zmm13, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm1, %zmm13, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm2{cdab}, %zmm13, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm2, %zmm13, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm3{cdab}, %zmm16, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm3, %zmm16, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm4{cdab}, %zmm16, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm4, %zmm16, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm5{cdab}, %zmm16, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm5, %zmm16, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm6{cdab}, %zmm19, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm6, %zmm19, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm7{cdab}, %zmm19, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm7, %zmm19, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm8{cdab}, %zmm19, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm8, %zmm19, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm9{cdab}, %zmm22, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm9, %zmm22, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm10{cdab}, %zmm22, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm10, %zmm22, %zmm29
	vmovnrngoapd  %zmm24,   64(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vmovapd        128(%rdx), %zmm24
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm11{cdab}, %zmm22, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm11, %zmm22, %zmm29
	vmovnrngoapd  %zmm25,  256(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vmovapd        320(%rdx), %zmm25
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm0{cdab}, %zmm14, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm0, %zmm14, %zmm29
	vmovnrngoapd  %zmm26,  448(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vmovapd        512(%rdx), %zmm26
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm1{cdab}, %zmm14, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm1, %zmm14, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm2{cdab}, %zmm14, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm2, %zmm14, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm3{cdab}, %zmm17, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm3, %zmm17, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm4{cdab}, %zmm17, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm4, %zmm17, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm5{cdab}, %zmm17, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm5, %zmm17, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm6{cdab}, %zmm20, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm6, %zmm20, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm7{cdab}, %zmm20, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm7, %zmm20, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm8{cdab}, %zmm20, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm8, %zmm20, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm9{cdab}, %zmm23, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm9, %zmm23, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm10{cdab}, %zmm23, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm10, %zmm23, %zmm29
	vmovnrngoapd  %zmm24,  128(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm11{cdab}, %zmm23, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm11, %zmm23, %zmm29
	vmovnrngoapd  %zmm25,  320(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vaddpd        %zmm26, %zmm29, %zmm26
	vmovnrngoapd  %zmm26,  512(%rdx)
	addq          $576, %rdx
	orq           %rcx, %rcx
	jle    _vmx_deriv_dag_lab5
	movq            24(%rbx), %rcx
	addq          %r12, %rcx
	addq          %r14, %rcx
	vmovapd          0(%rcx), %zmm0
	vmovapd         64(%rcx), %zmm1
	vmovapd        128(%rcx), %zmm2
	vmovapd        192(%rcx), %zmm3
	vmovapd        256(%rcx), %zmm4
	vmovapd        320(%rcx), %zmm5
	jmp    _vmx_deriv_dag_lab6
_vmx_deriv_dag_lab5:
	movq            24(%rbx), %rcx
	addq          %r8, %rcx
	addq          %r14, %rcx
	addq          %r14, %rcx
	vmovapd          0(%rcx), %zmm12
	vmovapd         64(%rcx), %zmm13
	vmovapd        128(%rcx), %zmm14
	vmovapd        192(%rcx), %zmm15
	vmovapd        256(%rcx), %zmm16
	vmovapd        320(%rcx), %zmm17
	vmovapd        384(%rcx), %zmm18
	vmovapd        448(%rcx), %zmm19
	vmovapd        512(%rcx), %zmm20
	vaddpd        %zmm18, %zmm15, %zmm3
	vmovapd        576(%rcx), %zmm21
	vaddpd        %zmm19, %zmm16, %zmm4
	vmovapd        640(%rcx), %zmm22
	vaddpd        %zmm20, %zmm17, %zmm5
	vmovapd        704(%rcx), %zmm23
	vsubpd        %zmm21, %zmm12, %zmm0
	vsubpd        %zmm22, %zmm13, %zmm1
	vsubpd        %zmm23, %zmm14, %zmm2
_vmx_deriv_dag_lab6:
	vaddpd        %zmm3, %zmm27, %zmm6
	vmovapd          0(%r9), %zmm12
	vsubpd        %zmm0, %zmm27, %zmm9
	vmovapd        192(%r9), %zmm15
	vaddpd        %zmm4, %zmm27, %zmm7
	vmovapd        384(%r9), %zmm18
	vsubpd        %zmm1, %zmm27, %zmm10
	vmovapd        576(%r9), %zmm21
	vaddpd        %zmm5, %zmm27, %zmm8
	vmovapd         64(%r9), %zmm13
	vsubpd        %zmm2, %zmm27, %zmm11
	vmovapd        256(%r9), %zmm16
	vmulpd        %zmm0, %zmm12, %zmm29
	vmovapd        448(%r9), %zmm19
	vmulpd        %zmm0{cdab}, %zmm12, %zmm30
	vmovapd        640(%r9), %zmm22
	vmovapd        128(%r9), %zmm14
	vmovapd        320(%r9), %zmm17
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vmovapd        512(%r9), %zmm20
	vmovapd        704(%r9), %zmm23
	vmovapd          0(%rdx), %zmm24
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmovapd        192(%rdx), %zmm25
	vmulpd        %zmm1{cdab}, %zmm12, %zmm30
	vmovapd        384(%rdx), %zmm26
	movq          %r15, %rcx
	andq          $32, %rcx
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm1, %zmm12, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm2{cdab}, %zmm12, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm2, %zmm12, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm3{cdab}, %zmm15, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm3, %zmm15, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm4{cdab}, %zmm15, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm4, %zmm15, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm5{cdab}, %zmm15, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm5, %zmm15, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm6{cdab}, %zmm18, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm6, %zmm18, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm7{cdab}, %zmm18, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm7, %zmm18, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm8{cdab}, %zmm18, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm8, %zmm18, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm9{cdab}, %zmm21, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm9, %zmm21, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm10{cdab}, %zmm21, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm10, %zmm21, %zmm29
	vmovnrngoapd  %zmm24,    0(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vmovapd         64(%rdx), %zmm24
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm11{cdab}, %zmm21, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm11, %zmm21, %zmm29
	vmovnrngoapd  %zmm25,  192(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vmovapd        256(%rdx), %zmm25
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm0{cdab}, %zmm13, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm0, %zmm13, %zmm29
	vmovnrngoapd  %zmm26,  384(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vmovapd        448(%rdx), %zmm26
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm1{cdab}, %zmm13, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm1, %zmm13, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm2{cdab}, %zmm13, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm2, %zmm13, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm3{cdab}, %zmm16, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm3, %zmm16, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm4{cdab}, %zmm16, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm4, %zmm16, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm5{cdab}, %zmm16, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm5, %zmm16, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm6{cdab}, %zmm19, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm6, %zmm19, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm7{cdab}, %zmm19, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm7, %zmm19, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm8{cdab}, %zmm19, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm8, %zmm19, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm9{cdab}, %zmm22, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm9, %zmm22, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm10{cdab}, %zmm22, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm10, %zmm22, %zmm29
	vmovnrngoapd  %zmm24,   64(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vmovapd        128(%rdx), %zmm24
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm11{cdab}, %zmm22, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm11, %zmm22, %zmm29
	vmovnrngoapd  %zmm25,  256(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vmovapd        320(%rdx), %zmm25
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm0{cdab}, %zmm14, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm0, %zmm14, %zmm29
	vmovnrngoapd  %zmm26,  448(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vmovapd        512(%rdx), %zmm26
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm1{cdab}, %zmm14, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm1, %zmm14, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm2{cdab}, %zmm14, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm2, %zmm14, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm3{cdab}, %zmm17, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm3, %zmm17, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm4{cdab}, %zmm17, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm4, %zmm17, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm5{cdab}, %zmm17, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm5, %zmm17, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm6{cdab}, %zmm20, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm6, %zmm20, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm7{cdab}, %zmm20, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm7, %zmm20, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm8{cdab}, %zmm20, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm8, %zmm20, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm9{cdab}, %zmm23, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm9, %zmm23, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm10{cdab}, %zmm23, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm10, %zmm23, %zmm29
	vmovnrngoapd  %zmm24,  128(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm11{cdab}, %zmm23, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm11, %zmm23, %zmm29
	vmovnrngoapd  %zmm25,  320(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vaddpd        %zmm26, %zmm29, %zmm26
	vmovnrngoapd  %zmm26,  512(%rdx)
	addq          $576, %rdx
	orq           %rcx, %rcx
	jle    _vmx_deriv_dag_lab7
	movq            40(%rbx), %rcx
	addq          %r12, %rcx
	addq          %r14, %rcx
	vmovapd          0(%rcx), %zmm0
	vmovapd         64(%rcx), %zmm1
	vmovapd        128(%rcx), %zmm2
	vmovapd        192(%rcx), %zmm3
	vmovapd        256(%rcx), %zmm4
	vmovapd        320(%rcx), %zmm5
	jmp    _vmx_deriv_dag_lab8
_vmx_deriv_dag_lab7:
	movq            40(%rbx), %rcx
	addq          %r8, %rcx
	addq          %r14, %rcx
	addq          %r14, %rcx
	vmovapd          0(%rcx), %zmm12
	vmovapd         64(%rcx), %zmm13
	vmovapd        128(%rcx), %zmm14
	vmovapd        192(%rcx), %zmm15
	vmovapd        256(%rcx), %zmm16
	vmovapd        320(%rcx), %zmm17
	vmovapd        384(%rcx), %zmm18
	vmovapd        448(%rcx), %zmm19
	vmovapd        512(%rcx), %zmm20
	vsubpd        %zmm18{cdab}, %zmm12, %zmm0{%k7}
	vmovapd        576(%rcx), %zmm21
	vaddpd        %zmm18{cdab}, %zmm12, %zmm0{%k6}
	vmovapd        640(%rcx), %zmm22
	vsubpd        %zmm19{cdab}, %zmm13, %zmm1{%k7}
	vmovapd        704(%rcx), %zmm23
	vaddpd        %zmm19{cdab}, %zmm13, %zmm1{%k6}
	vsubpd        %zmm20{cdab}, %zmm14, %zmm2{%k7}
	vaddpd        %zmm20{cdab}, %zmm14, %zmm2{%k6}
	vaddpd        %zmm21{cdab}, %zmm15, %zmm3{%k7}
	vsubpd        %zmm21{cdab}, %zmm15, %zmm3{%k6}
	vaddpd        %zmm22{cdab}, %zmm16, %zmm4{%k7}
	vsubpd        %zmm22{cdab}, %zmm16, %zmm4{%k6}
	vaddpd        %zmm23{cdab}, %zmm17, %zmm5{%k7}
	vsubpd        %zmm23{cdab}, %zmm17, %zmm5{%k6}
_vmx_deriv_dag_lab8:
	vaddpd        %zmm0{cdab}, %zmm27, %zmm6{%k7}
	vmovapd          0(%r9), %zmm12
	vsubpd        %zmm0{cdab}, %zmm27, %zmm6{%k6}
	vmovapd        192(%r9), %zmm15
	vsubpd        %zmm3{cdab}, %zmm27, %zmm9{%k7}
	vmovapd        384(%r9), %zmm18
	vaddpd        %zmm3{cdab}, %zmm27, %zmm9{%k6}
	vmovapd        576(%r9), %zmm21
	vaddpd        %zmm1{cdab}, %zmm27, %zmm7{%k7}
	vmovapd         64(%r9), %zmm13
	vsubpd        %zmm1{cdab}, %zmm27, %zmm7{%k6}
	vmovapd        256(%r9), %zmm16
	vsubpd        %zmm4{cdab}, %zmm27, %zmm10{%k7}
	vmovapd        448(%r9), %zmm19
	vaddpd        %zmm4{cdab}, %zmm27, %zmm10{%k6}
	vmovapd        640(%r9), %zmm22
	vaddpd        %zmm2{cdab}, %zmm27, %zmm8{%k7}
	vmovapd        128(%r9), %zmm14
	vsubpd        %zmm2{cdab}, %zmm27, %zmm8{%k6}
	vmovapd        320(%r9), %zmm17
	vsubpd        %zmm5{cdab}, %zmm27, %zmm11{%k7}
	vmovapd        512(%r9), %zmm20
	vaddpd        %zmm5{cdab}, %zmm27, %zmm11{%k6}
	vmovapd        704(%r9), %zmm23
	vmulpd        %zmm0, %zmm12, %zmm29
	vmovapd          0(%rdx), %zmm24
	vmulpd        %zmm0{cdab}, %zmm12, %zmm30
	vmovapd        192(%rdx), %zmm25
	vmovapd        384(%rdx), %zmm26
	movq          %r15, %rcx
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	andq          $128, %rcx
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm1{cdab}, %zmm12, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm1, %zmm12, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm2{cdab}, %zmm12, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm2, %zmm12, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm3{cdab}, %zmm15, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm3, %zmm15, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm4{cdab}, %zmm15, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm4, %zmm15, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm5{cdab}, %zmm15, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm5, %zmm15, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm6{cdab}, %zmm18, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm6, %zmm18, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm7{cdab}, %zmm18, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm7, %zmm18, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm8{cdab}, %zmm18, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm8, %zmm18, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm9{cdab}, %zmm21, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm9, %zmm21, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm10{cdab}, %zmm21, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm10, %zmm21, %zmm29
	vmovnrngoapd  %zmm24,    0(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vmovapd         64(%rdx), %zmm24
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm11{cdab}, %zmm21, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm11, %zmm21, %zmm29
	vmovnrngoapd  %zmm25,  192(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vmovapd        256(%rdx), %zmm25
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm0{cdab}, %zmm13, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm0, %zmm13, %zmm29
	vmovnrngoapd  %zmm26,  384(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vmovapd        448(%rdx), %zmm26
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm1{cdab}, %zmm13, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm1, %zmm13, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm2{cdab}, %zmm13, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm2, %zmm13, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm3{cdab}, %zmm16, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm3, %zmm16, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm4{cdab}, %zmm16, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm4, %zmm16, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm5{cdab}, %zmm16, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm5, %zmm16, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm6{cdab}, %zmm19, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm6, %zmm19, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm7{cdab}, %zmm19, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm7, %zmm19, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm8{cdab}, %zmm19, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm8, %zmm19, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm9{cdab}, %zmm22, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm9, %zmm22, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm10{cdab}, %zmm22, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm10, %zmm22, %zmm29
	vmovnrngoapd  %zmm24,   64(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vmovapd        128(%rdx), %zmm24
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm11{cdab}, %zmm22, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm11, %zmm22, %zmm29
	vmovnrngoapd  %zmm25,  256(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vmovapd        320(%rdx), %zmm25
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm0{cdab}, %zmm14, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm0, %zmm14, %zmm29
	vmovnrngoapd  %zmm26,  448(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vmovapd        512(%rdx), %zmm26
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm1{cdab}, %zmm14, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm1, %zmm14, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm2{cdab}, %zmm14, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm2, %zmm14, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm3{cdab}, %zmm17, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm3, %zmm17, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm4{cdab}, %zmm17, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm4, %zmm17, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm5{cdab}, %zmm17, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm5, %zmm17, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm6{cdab}, %zmm20, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm6, %zmm20, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm7{cdab}, %zmm20, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm7, %zmm20, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm8{cdab}, %zmm20, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm8, %zmm20, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm9{cdab}, %zmm23, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm9, %zmm23, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm10{cdab}, %zmm23, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm10, %zmm23, %zmm29
	vmovnrngoapd  %zmm24,  128(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm11{cdab}, %zmm23, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm11, %zmm23, %zmm29
	vmovnrngoapd  %zmm25,  320(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vaddpd        %zmm26, %zmm29, %zmm26
	vmovnrngoapd  %zmm26,  512(%rdx)
	addq          $576, %rdx
	orq           %rcx, %rcx
	jle    _vmx_deriv_dag_lab9
	movq            56(%rbx), %rcx
	addq          %r12, %rcx
	addq          %r14, %rcx
	vmovapd          0(%rcx), %zmm0
	vmovapd         64(%rcx), %zmm1
	vmovapd        128(%rcx), %zmm2
	vmovapd        192(%rcx), %zmm3
	vmovapd        256(%rcx), %zmm4
	vmovapd        320(%rcx), %zmm5
	jmp    _vmx_deriv_dag_lab10
_vmx_deriv_dag_lab9:
	movq            56(%rbx), %rcx
	addq          %r8, %rcx
	addq          %r14, %rcx
	addq          %r14, %rcx
	vmovapd          0(%rcx), %zmm12
	vmovapd         64(%rcx), %zmm13
	vmovapd        128(%rcx), %zmm14
	vmovapd        192(%rcx), %zmm15
	vmovapd        256(%rcx), %zmm16
	vmovapd        320(%rcx), %zmm17
	vmovapd        384(%rcx), %zmm18
	vmovapd        448(%rcx), %zmm19
	vmovapd        512(%rcx), %zmm20
	vaddpd        %zmm18, %zmm12, %zmm0
	vmovapd        576(%rcx), %zmm21
	vaddpd        %zmm19, %zmm13, %zmm1
	vmovapd        640(%rcx), %zmm22
	vaddpd        %zmm20, %zmm14, %zmm2
	vmovapd        704(%rcx), %zmm23
	vaddpd        %zmm21, %zmm15, %zmm3
	vaddpd        %zmm22, %zmm16, %zmm4
	vaddpd        %zmm23, %zmm17, %zmm5
_vmx_deriv_dag_lab10:
	vaddpd        %zmm0, %zmm27, %zmm6
	vmovapd          0(%r9), %zmm12
	vaddpd        %zmm3, %zmm27, %zmm9
	vmovapd        192(%r9), %zmm15
	vaddpd        %zmm1, %zmm27, %zmm7
	vmovapd        384(%r9), %zmm18
	vaddpd        %zmm4, %zmm27, %zmm10
	vmovapd        576(%r9), %zmm21
	vaddpd        %zmm2, %zmm27, %zmm8
	vmovapd         64(%r9), %zmm13
	vaddpd        %zmm5, %zmm27, %zmm11
	vmovapd        256(%r9), %zmm16
	vmulpd        %zmm0, %zmm12, %zmm29
	vmovapd        448(%r9), %zmm19
	vmulpd        %zmm0{cdab}, %zmm12, %zmm30
	vmovapd        640(%r9), %zmm22
	vmovapd        128(%r9), %zmm14
	vmovapd        320(%r9), %zmm17
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vmovapd        512(%r9), %zmm20
	vmovapd        704(%r9), %zmm23
	vmovapd          0(%rdx), %zmm24
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmovapd        192(%rdx), %zmm25
	vmulpd        %zmm1{cdab}, %zmm12, %zmm30
	vmovapd        384(%rdx), %zmm26
	addq          $768, %r9
	addq          $384, %r14
	addq          $-1, %r11
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm1, %zmm12, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm2{cdab}, %zmm12, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm2, %zmm12, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm3{cdab}, %zmm15, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm3, %zmm15, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm4{cdab}, %zmm15, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm4, %zmm15, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm5{cdab}, %zmm15, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm5, %zmm15, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm6{cdab}, %zmm18, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm6, %zmm18, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm7{cdab}, %zmm18, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm7, %zmm18, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm8{cdab}, %zmm18, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm8, %zmm18, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm9{cdab}, %zmm21, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm9, %zmm21, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm10{cdab}, %zmm21, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm10, %zmm21, %zmm29
	vmovnrngoapd  %zmm24,    0(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vmovapd         64(%rdx), %zmm24
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm11{cdab}, %zmm21, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm11, %zmm21, %zmm29
	vmovnrngoapd  %zmm25,  192(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vmovapd        256(%rdx), %zmm25
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm0{cdab}, %zmm13, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm0, %zmm13, %zmm29
	vmovnrngoapd  %zmm26,  384(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vmovapd        448(%rdx), %zmm26
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm1{cdab}, %zmm13, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm1, %zmm13, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm2{cdab}, %zmm13, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm2, %zmm13, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm3{cdab}, %zmm16, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm3, %zmm16, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm4{cdab}, %zmm16, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm4, %zmm16, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm5{cdab}, %zmm16, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm5, %zmm16, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm6{cdab}, %zmm19, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm6, %zmm19, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm7{cdab}, %zmm19, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm7, %zmm19, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm8{cdab}, %zmm19, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm8, %zmm19, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm9{cdab}, %zmm22, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm9, %zmm22, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm10{cdab}, %zmm22, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm10, %zmm22, %zmm29
	vmovnrngoapd  %zmm24,   64(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vmovapd        128(%rdx), %zmm24
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm11{cdab}, %zmm22, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm11, %zmm22, %zmm29
	vmovnrngoapd  %zmm25,  256(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vmovapd        320(%rdx), %zmm25
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm0{cdab}, %zmm14, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm0, %zmm14, %zmm29
	vmovnrngoapd  %zmm26,  448(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vmovapd        512(%rdx), %zmm26
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm1{cdab}, %zmm14, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm1, %zmm14, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm2{cdab}, %zmm14, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm2, %zmm14, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm3{cdab}, %zmm17, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm3, %zmm17, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm4{cdab}, %zmm17, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm4, %zmm17, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm5{cdab}, %zmm17, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm5, %zmm17, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm6{cdab}, %zmm20, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm6, %zmm20, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm7{cdab}, %zmm20, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm7, %zmm20, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm8{cdab}, %zmm20, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm8, %zmm20, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm9{cdab}, %zmm23, %zmm30
	vaddpd        %zmm26, %zmm29, %zmm26
	vmulpd        %zmm9, %zmm23, %zmm29
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm10{cdab}, %zmm23, %zmm30
	vaddpd        %zmm24, %zmm29, %zmm24
	vmulpd        %zmm10, %zmm23, %zmm29
	vmovnrngoapd  %zmm24,  128(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vmulpd        %zmm11{cdab}, %zmm23, %zmm30
	vaddpd        %zmm25, %zmm29, %zmm25
	vmulpd        %zmm11, %zmm23, %zmm29
	vmovnrngoapd  %zmm25,  320(%rdx)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm29{%k7}
	vsubrpd       %zmm30{cdab}, %zmm30, %zmm29{%k6}
	vaddpd        %zmm26, %zmm29, %zmm26
	vmovnrngoapd  %zmm26,  512(%rdx)
	addq          $576, %rdx
	orq           %r11, %r11
	jg     _vmx_deriv_dag_lab2
	movq          %rdx, %rsi
	addq          $128, %rbx
	addq          $-1, %rax
	orq           %rax, %rax
	jg     _vmx_deriv_dag_lab1
_vmx_deriv_dag_lab0:
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
.type vmx_deriv_dag,@function
.size vmx_deriv_dag,.-vmx_deriv_dag
.data
