/* 
 * BAGEL machine generated output.  
 *   
 * It is provided under the GNU pubic License V2  
 * It is provided as is, and is not guaranteed fit for any purpose.
 * BAGEL was written by Peter Boyle  
 */  
.text
.align 16,0x90
.globl vmx_merge
vmx_merge:
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
	movq          %rdx, %r10
	movq          %rbx, 2176(%rsp)
	movq          %rbp, 2184(%rsp)
	movq          %r12, 2208(%rsp)
	movq          %r13, 2216(%rsp)
	movq          %r14, 2224(%rsp)
	movq          %r15, 2232(%rsp)
	addq          $0, %rsp
	movq          %rdi, %rbx
	movq          %rsi, %rbp
	movq          %rcx, %r11
_vmx_merge_lab0:
	orq           %r11, %r11
	jle    _vmx_merge_lab1
_vmx_merge_lab2:
	vmovapd          0(%rbp), %zmm0
	vmovapd          0(%r10), %zmm3
	vmovapd         64(%rbp), %zmm1
	vpermf32x4    $78, %zmm0, %zmm7
	vmovapd         64(%r10), %zmm4
	vpermf32x4    $78, %zmm3, %zmm6
	vmovapd        128(%rbp), %zmm2
	vpermf32x4    $78, %zmm1, %zmm9
	vmovapd        128(%r10), %zmm5
	vpermf32x4    $78, %zmm4, %zmm8
	addq          $192, %rbp
	vpermf32x4    $78, %zmm2, %zmm11
	addq          $192, %r10
	vpermf32x4    $78, %zmm5, %zmm10
	addq          $-1, %r11
	vblendmpd     %zmm3, %zmm7, %zmm7{%k2}
	vblendmpd     %zmm0, %zmm6, %zmm6{%k2}
	vblendmpd     %zmm4, %zmm9, %zmm9{%k2}
	vblendmpd     %zmm1, %zmm8, %zmm8{%k2}
	vblendmpd     %zmm5, %zmm11, %zmm11{%k2}
	vblendmpd     %zmm2, %zmm10, %zmm10{%k2}
	vmovnrngoapd  %zmm7,   64(%rbx)
	vmovnrngoapd  %zmm6,    0(%rbx)
	vmovnrngoapd  %zmm9,  192(%rbx)
	vmovnrngoapd  %zmm8,  128(%rbx)
	vmovnrngoapd  %zmm11,  320(%rbx)
	vmovnrngoapd  %zmm10,  256(%rbx)
	addq          $384, %rbx
	orq           %r11, %r11
	jg     _vmx_merge_lab2
_vmx_merge_lab1:
	addq          $0, %rsp
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
.type vmx_merge,@function
.size vmx_merge,.-vmx_merge
.data
