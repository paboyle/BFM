/* 
 * BAGEL machine generated output.  
 *   
 * It is provided under the GNU pubic License V2  
 * It is provided as is, and is not guaranteed fit for any purpose.
 * BAGEL was written by Peter Boyle  
 */  
.text
.align 16,0x90
.globl vmx_cg_s
vmx_cg_s:
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
	jle    _vmx_cg_s_lab0
_vmx_cg_s_lab1:
	vmovaps          0(%rcx), %zmm2
	vmovaps          0(%r8), %zmm8
	vmovaps          0(%r9), %zmm14
	vmovaps          0(%rax), %zmm20
	vmovaps         64(%rcx), %zmm4
	vmovaps         64(%r8), %zmm10
	vmovaps         64(%r9), %zmm16
	vmovaps         64(%rax), %zmm22
	vmovaps        128(%rcx), %zmm6
	vmovaps        128(%r8), %zmm12
	vmovaps        128(%r9), %zmm18
	vmovaps        128(%rax), %zmm24
	movq          %rcx, %r11
	movq          %r8, %r12
	movq          %r9, %r13
	addq          $192, %rcx
	addq          $192, %r8
	addq          $192, %r9
	vmovaps         32(%r11), %zmm3
	vmovaps         32(%r12), %zmm9
	vmovaps         32(%r13), %zmm15
	vmovaps         32(%rax), %zmm21
	vmovaps         96(%r11), %zmm5
	vmovaps         96(%r12), %zmm11
	vmovaps         96(%r13), %zmm17
	vmovaps         96(%rax), %zmm23
	vmovaps        160(%r11), %zmm7
	vmovaps        160(%r12), %zmm13
	vmovaps        160(%r13), %zmm19
	vmovaps        160(%rax), %zmm25
	vfnmadd231pd  %zmm20, %zmm0, %zmm14
	vfnmadd231pd  %zmm21, %zmm0, %zmm15
	vfnmadd231pd  %zmm22, %zmm0, %zmm16
	vfnmadd231pd  %zmm23, %zmm0, %zmm17
	vfnmadd231pd  %zmm24, %zmm0, %zmm18
	vfnmadd231pd  %zmm25, %zmm0, %zmm19
	vfmadd231pd   %zmm8, %zmm0, %zmm2
	vmovnrngoaps  %zmm14,    0(%r13)
	vfmadd231pd   %zmm9, %zmm0, %zmm3
	vmovnrngoaps  %zmm15,   32(%r13)
	vfmadd231pd   %zmm10, %zmm0, %zmm4
	vmovnrngoaps  %zmm16,   64(%r13)
	vfmadd231pd   %zmm11, %zmm0, %zmm5
	vmovnrngoaps  %zmm17,   96(%r13)
	vfmadd231pd   %zmm12, %zmm0, %zmm6
	vmovnrngoaps  %zmm18,  128(%r13)
	vfmadd231pd   %zmm13, %zmm0, %zmm7
	vmovnrngoaps  %zmm19,  160(%r13)
	addq          $192, %rax
	vfmadd213pd   %zmm14, %zmm1, %zmm8
	vmovnrngoaps  %zmm2,    0(%r11)
	vfmadd213pd   %zmm15, %zmm1, %zmm9
	vmovnrngoaps  %zmm3,   32(%r11)
	vfmadd213pd   %zmm16, %zmm1, %zmm10
	vmovnrngoaps  %zmm4,   64(%r11)
	vfmadd213pd   %zmm17, %zmm1, %zmm11
	vmovnrngoaps  %zmm5,   96(%r11)
	vfmadd213pd   %zmm18, %zmm1, %zmm12
	vmovnrngoaps  %zmm6,  128(%r11)
	vfmadd213pd   %zmm19, %zmm1, %zmm13
	vmovnrngoaps  %zmm7,  160(%r11)
	vfmadd231pd   %zmm14, %zmm14, %zmm26
	vmovnrngoaps  %zmm8,    0(%r12)
	vfmadd231pd   %zmm15, %zmm15, %zmm26
	vmovnrngoaps  %zmm9,   32(%r12)
	vfmadd231pd   %zmm16, %zmm16, %zmm26
	vmovnrngoaps  %zmm10,   64(%r12)
	vfmadd231pd   %zmm17, %zmm17, %zmm26
	vmovnrngoaps  %zmm11,   96(%r12)
	vfmadd231pd   %zmm18, %zmm18, %zmm26
	vmovnrngoaps  %zmm12,  128(%r12)
	vfmadd231pd   %zmm19, %zmm19, %zmm26
	vmovnrngoaps  %zmm13,  160(%r12)
	addq          $-1, %rbx
	orq           %rbx, %rbx
	jg     _vmx_cg_s_lab1
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
_vmx_cg_s_lab0:
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
.type vmx_cg_s,@function
.size vmx_cg_s,.-vmx_cg_s
.data
