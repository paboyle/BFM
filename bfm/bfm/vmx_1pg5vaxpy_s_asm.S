/* 
 * BAGEL machine generated output.  
 *   
 * It is provided under the GNU pubic License V2  
 * It is provided as is, and is not guaranteed fit for any purpose.
 * BAGEL was written by Peter Boyle  
 */  
.text
.align 16,0x90
.globl vmx_1pG5vaxpy_s
vmx_1pG5vaxpy_s:
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
	movq          %rdi, %rbp
	movq          %rsi, %r13
	movq          %rdx, %rbx
	movq          %rcx, %r12
	vbroadcastsd     0(%r13), %zmm0
	orq           %r12, %r12
	jle    _vmx_1pG5vaxpy_s_lab0
	vbroadcastss     0(%rbx), %zmm1
	vbroadcastss     0(%rbp), %zmm13
	vbroadcastss     4(%rbx), %zmm2
	vbroadcastss     4(%rbp), %zmm14
	vbroadcastss     8(%rbx), %zmm3
	vbroadcastss     8(%rbp), %zmm15
	vbroadcastss    12(%rbx), %zmm4
	vbroadcastss    12(%rbp), %zmm16
	vbroadcastss    16(%rbx), %zmm5
	vbroadcastss    16(%rbp), %zmm17
	vbroadcastss    20(%rbx), %zmm6
	vbroadcastss    20(%rbp), %zmm18
_vmx_1pG5vaxpy_s_lab1:
	vmovapd       %zmm0, %zmm25{%k5}
	vbroadcastss    24(%rbx), %zmm7
	vfmadd213pd   %zmm13, %zmm1, %zmm25{%k5}
	vbroadcastss    24(%rbp), %zmm19
	vmovapd       %zmm0, %zmm26{%k5}
	vbroadcastss    28(%rbx), %zmm8
	vfmadd213pd   %zmm14, %zmm2, %zmm26{%k5}
	vbroadcastss    28(%rbp), %zmm20
	vmovapd       %zmm0, %zmm27{%k5}
	vbroadcastss    32(%rbx), %zmm9
	vfmadd213pd   %zmm15, %zmm3, %zmm27{%k5}
	vbroadcastss    32(%rbp), %zmm21
	vmovapd       %zmm0, %zmm28{%k5}
	vbroadcastss    36(%rbx), %zmm10
	vfmadd213pd   %zmm16, %zmm4, %zmm28{%k5}
	vbroadcastss    36(%rbp), %zmm22
	vmovapd       %zmm0, %zmm29{%k5}
	vbroadcastss    40(%rbx), %zmm11
	vfmadd213pd   %zmm17, %zmm5, %zmm29{%k5}
	vbroadcastss    40(%rbp), %zmm23
	vmovapd       %zmm0, %zmm30{%k5}
	vbroadcastss    44(%rbx), %zmm12
	vfmadd213pd   %zmm18, %zmm6, %zmm30{%k5}
	vbroadcastss    44(%rbp), %zmm24
	vpackstorelps %zmm25,    0(%rbp){%k5}
	vmovapd       %zmm0, %zmm25{%k5}
	vpackstorelps %zmm26,    4(%rbp){%k5}
	vfmadd213pd   %zmm19, %zmm7, %zmm25{%k5}
	addq          $96, %rbx
	vmovapd       %zmm0, %zmm26{%k5}
	vpackstorelps %zmm27,    8(%rbp){%k5}
	vfmadd213pd   %zmm20, %zmm8, %zmm26{%k5}
	vpackstorelps %zmm28,   12(%rbp){%k5}
	vmovapd       %zmm0, %zmm27{%k5}
	vbroadcastss     0(%rbx), %zmm1
	vfmadd213pd   %zmm21, %zmm9, %zmm27{%k5}
	vpackstorelps %zmm29,   16(%rbp){%k5}
	vmovapd       %zmm0, %zmm28{%k5}
	vpackstorelps %zmm30,   20(%rbp){%k5}
	vfmadd213pd   %zmm22, %zmm10, %zmm28{%k5}
	vbroadcastss     4(%rbx), %zmm2
	vmovapd       %zmm0, %zmm29{%k5}
	vpackstorelps %zmm25,   24(%rbp){%k5}
	vfmadd213pd   %zmm23, %zmm11, %zmm29{%k5}
	vpackstorelps %zmm26,   28(%rbp){%k5}
	vmovapd       %zmm0, %zmm30{%k5}
	vbroadcastss     8(%rbx), %zmm3
	vfmadd213pd   %zmm24, %zmm12, %zmm30{%k5}
	vpackstorelps %zmm27,   32(%rbp){%k5}
	vpackstorelps %zmm28,   36(%rbp){%k5}
	vbroadcastss    12(%rbx), %zmm4
	vpackstorelps %zmm29,   40(%rbp){%k5}
	vpackstorelps %zmm30,   44(%rbp){%k5}
	addq          $96, %rbp
	vbroadcastss    16(%rbx), %zmm5
	vbroadcastss     0(%rbp), %zmm13
	vbroadcastss     4(%rbp), %zmm14
	vbroadcastss     8(%rbp), %zmm15
	vbroadcastss    12(%rbp), %zmm16
	vbroadcastss    16(%rbp), %zmm17
	vbroadcastss    20(%rbx), %zmm6
	vbroadcastss    20(%rbp), %zmm18
	addq          $-1, %r12
	orq           %r12, %r12
	jg     _vmx_1pG5vaxpy_s_lab1
_vmx_1pG5vaxpy_s_lab0:
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
.type vmx_1pG5vaxpy_s,@function
.size vmx_1pG5vaxpy_s,.-vmx_1pG5vaxpy_s
.data
