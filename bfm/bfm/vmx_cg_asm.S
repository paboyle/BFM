/* 
 * BAGEL machine generated output.  
 *   
 * It is provided under the GNU pubic License V2  
 * It is provided as is, and is not guaranteed fit for any purpose.
 * BAGEL was written by Peter Boyle  
 */  
.text
.align 16,0x90
.globl vmx_cg
vmx_cg:
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

	addq          $-2240, %rsp
	movq          %rbx, 2176(%rsp)
	movq          %rbp, 2184(%rsp)
	movq          %r12, 2208(%rsp)
	movq          %r13, 2216(%rsp)
	movq          %r14, 2224(%rsp)
	movq          %r15, 2232(%rsp)
	movq          %rdi, %r15
	movq          %rsp, %rdi
	movq             0(%r15), %rbp
	addq          $8, %r15
	movq             0(%r15), %r10
	addq          $8, %r15
	movq             0(%r15), %rcx
	addq          $8, %r15
	movq             0(%r15), %r8
	addq          $8, %r15
	movq             0(%r15), %r9
	addq          $8, %r15
	movq             0(%r15), %rax
	addq          $8, %r15
	movq             0(%r15), %rbx
	addq          $8, %r15
	movq             0(%r15), %r14
	addq          $8, %r15
	vbroadcastsd     0(%r14), %zmm26
	vbroadcastsd     0(%rbp), %zmm0
	vbroadcastsd     0(%r10), %zmm1
	orq           %rbx, %rbx
	jle    _vmx_cg_lab0
_vmx_cg_lab1:
	vmovapd          0(%rcx), %zmm2
	vmovapd          0(%r8), %zmm8
	vmovapd          0(%r9), %zmm14
	vmovapd          0(%rax), %zmm20
	vmovapd        128(%rcx), %zmm4
	vmovapd        128(%r8), %zmm10
	vmovapd        128(%r9), %zmm16
	vmovapd        128(%rax), %zmm22
	vmovapd        256(%rcx), %zmm6
	vmovapd        256(%r8), %zmm12
	vmovapd        256(%r9), %zmm18
	vmovapd        256(%rax), %zmm24
	movq          %rcx, %r11
	movq          %r8, %r12
	movq          %r9, %r13
	addq          $384, %rcx
	addq          $384, %r8
	addq          $384, %r9
	vmovapd         64(%r11), %zmm3
	vmovapd         64(%r12), %zmm9
	vmovapd         64(%r13), %zmm15
	vmovapd         64(%rax), %zmm21
	vmovapd        192(%r11), %zmm5
	vmovapd        192(%r12), %zmm11
	vmovapd        192(%r13), %zmm17
	vmovapd        192(%rax), %zmm23
	vmovapd        320(%r11), %zmm7
	vmovapd        320(%r12), %zmm13
	vmovapd        320(%r13), %zmm19
	vmovapd        320(%rax), %zmm25
	vfnmadd231pd  %zmm20, %zmm0, %zmm14
	vfnmadd231pd  %zmm21, %zmm0, %zmm15
	vfnmadd231pd  %zmm22, %zmm0, %zmm16
	vfnmadd231pd  %zmm23, %zmm0, %zmm17
	vfnmadd231pd  %zmm24, %zmm0, %zmm18
	vfnmadd231pd  %zmm25, %zmm0, %zmm19
	vfmadd231pd   %zmm8, %zmm0, %zmm2
	vmovnrngoapd  %zmm14,    0(%r13)
	vfmadd231pd   %zmm9, %zmm0, %zmm3
	vmovnrngoapd  %zmm15,   64(%r13)
	vfmadd231pd   %zmm10, %zmm0, %zmm4
	vmovnrngoapd  %zmm16,  128(%r13)
	vfmadd231pd   %zmm11, %zmm0, %zmm5
	vmovnrngoapd  %zmm17,  192(%r13)
	vfmadd231pd   %zmm12, %zmm0, %zmm6
	vmovnrngoapd  %zmm18,  256(%r13)
	vfmadd231pd   %zmm13, %zmm0, %zmm7
	vmovnrngoapd  %zmm19,  320(%r13)
	addq          $384, %rax
	vfmadd213pd   %zmm14, %zmm1, %zmm8
	vmovnrngoapd  %zmm2,    0(%r11)
	vfmadd213pd   %zmm15, %zmm1, %zmm9
	vmovnrngoapd  %zmm3,   64(%r11)
	vfmadd213pd   %zmm16, %zmm1, %zmm10
	vmovnrngoapd  %zmm4,  128(%r11)
	vfmadd213pd   %zmm17, %zmm1, %zmm11
	vmovnrngoapd  %zmm5,  192(%r11)
	vfmadd213pd   %zmm18, %zmm1, %zmm12
	vmovnrngoapd  %zmm6,  256(%r11)
	vfmadd213pd   %zmm19, %zmm1, %zmm13
	vmovnrngoapd  %zmm7,  320(%r11)
	vfmadd231pd   %zmm14, %zmm14, %zmm26
	vmovnrngoapd  %zmm8,    0(%r12)
	vfmadd231pd   %zmm15, %zmm15, %zmm26
	vmovnrngoapd  %zmm9,   64(%r12)
	vfmadd231pd   %zmm16, %zmm16, %zmm26
	vmovnrngoapd  %zmm10,  128(%r12)
	vfmadd231pd   %zmm17, %zmm17, %zmm26
	vmovnrngoapd  %zmm11,  192(%r12)
	vfmadd231pd   %zmm18, %zmm18, %zmm26
	vmovnrngoapd  %zmm12,  256(%r12)
	vfmadd231pd   %zmm19, %zmm19, %zmm26
	vmovnrngoapd  %zmm13,  320(%r12)
	addq          $-1, %rbx
	orq           %rbx, %rbx
	jg     _vmx_cg_lab1
	vmovnrngoapd  %zmm26,    0(%rdi)
	vbroadcastsd     0(%rdi), %zmm28
	vbroadcastsd     8(%rdi), %zmm27
	vaddpd        %zmm27, %zmm28, %zmm28{%k5}
	vbroadcastsd    16(%rdi), %zmm27
	vaddpd        %zmm27, %zmm28, %zmm28{%k5}
	vbroadcastsd    24(%rdi), %zmm27
	vaddpd        %zmm27, %zmm28, %zmm28{%k5}
	vbroadcastsd    32(%rdi), %zmm27
	vaddpd        %zmm27, %zmm28, %zmm28{%k5}
	vbroadcastsd    40(%rdi), %zmm27
	vaddpd        %zmm27, %zmm28, %zmm28{%k5}
	vbroadcastsd    48(%rdi), %zmm27
	vaddpd        %zmm27, %zmm28, %zmm28{%k5}
	vbroadcastsd    56(%rdi), %zmm27
	vaddpd        %zmm27, %zmm28, %zmm28{%k5}
	vpackstorelpd %zmm28,    0(%r14){%k5}
_vmx_cg_lab0:
	movq          2232(%rsp), %r15
	movq          2224(%rsp), %r14
	movq          2216(%rsp), %r13
	movq          2208(%rsp), %r12
	movq          2184(%rsp), %rbp
	movq          2176(%rsp), %rbx
	addq          $2240, %rsp

	# ----- Exit_Routine -----
	movq     (%rsp), %rsp
	popq     %r10
	ret
.align 16,0x90
.type vmx_cg,@function
.size vmx_cg,.-vmx_cg
.data
