/* 
 * BAGEL machine generated output.  
 *   
 * It is provided under the GNU pubic License V2  
 * It is provided as is, and is not guaranteed fit for any purpose.
 * BAGEL was written by Peter Boyle  
 */  
.text
.align 16,0x90
.globl vmx_tmass
vmx_tmass:
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
	movq          %rsi, %r11
	movq          %rdx, %r12
	movq          %rcx, %rbx
	movq          %r8, %r10
	vbroadcastsd     0(%r11), %zmm0
	vbroadcastsd     0(%r12), %zmm1
_vmx_tmass_lab0:
	orq           %r10, %r10
	jle    _vmx_tmass_lab1
_vmx_tmass_lab2:
	vmovapd          0(%rbx), %zmm2
	vmovapd        128(%rbx), %zmm4
	vmovapd         64(%rbx), %zmm3
	vmulpd        %zmm2, %zmm0, %zmm5
	vmulpd        %zmm2, %zmm1, %zmm2
	vsubpd        %zmm2{cdab}, %zmm5, %zmm5{%k7}
	vaddpd        %zmm2{cdab}, %zmm5, %zmm5{%k6}
	vmulpd        %zmm4, %zmm0, %zmm7
	vmulpd        %zmm4, %zmm1, %zmm4
	vsubpd        %zmm4{cdab}, %zmm7, %zmm7{%k7}
	vaddpd        %zmm4{cdab}, %zmm7, %zmm7{%k6}
	vmulpd        %zmm3, %zmm0, %zmm6
	vmulpd        %zmm3, %zmm1, %zmm3
	vsubpd        %zmm3{cdab}, %zmm6, %zmm6{%k7}
	vaddpd        %zmm3{cdab}, %zmm6, %zmm6{%k6}
	vmovnrngoapd  %zmm5,    0(%rbp)
	vmovnrngoapd  %zmm7,  128(%rbp)
	vmovnrngoapd  %zmm6,   64(%rbp)
	addq          $192, %rbx
	addq          $192, %rbp
	vmovapd          0(%rbx), %zmm2
	vmovapd        128(%rbx), %zmm4
	vmovapd         64(%rbx), %zmm3
	vmulpd        %zmm2, %zmm0, %zmm5
	vmulpd        %zmm2, %zmm1, %zmm2
	vsubpd        %zmm2{cdab}, %zmm5, %zmm5{%k7}
	vaddpd        %zmm2{cdab}, %zmm5, %zmm5{%k6}
	vmulpd        %zmm4, %zmm0, %zmm7
	vmulpd        %zmm4, %zmm1, %zmm4
	vsubpd        %zmm4{cdab}, %zmm7, %zmm7{%k7}
	vaddpd        %zmm4{cdab}, %zmm7, %zmm7{%k6}
	vmulpd        %zmm3, %zmm0, %zmm6
	vmulpd        %zmm3, %zmm1, %zmm3
	vsubpd        %zmm3{cdab}, %zmm6, %zmm6{%k7}
	vaddpd        %zmm3{cdab}, %zmm6, %zmm6{%k6}
	vmovnrngoapd  %zmm5,    0(%rbp)
	vmovnrngoapd  %zmm7,  128(%rbp)
	vmovnrngoapd  %zmm6,   64(%rbp)
	addq          $192, %rbx
	addq          $192, %rbp
	vmovapd          0(%rbx), %zmm2
	vmovapd        128(%rbx), %zmm4
	vmovapd         64(%rbx), %zmm3
	vmulpd        %zmm2, %zmm0, %zmm5
	vmulpd        %zmm2, %zmm1, %zmm2
	vaddpd        %zmm2{cdab}, %zmm5, %zmm5{%k7}
	vsubpd        %zmm2{cdab}, %zmm5, %zmm5{%k6}
	vmulpd        %zmm4, %zmm0, %zmm7
	vmulpd        %zmm4, %zmm1, %zmm4
	vaddpd        %zmm4{cdab}, %zmm7, %zmm7{%k7}
	vsubpd        %zmm4{cdab}, %zmm7, %zmm7{%k6}
	vmulpd        %zmm3, %zmm0, %zmm6
	vmulpd        %zmm3, %zmm1, %zmm3
	vaddpd        %zmm3{cdab}, %zmm6, %zmm6{%k7}
	vsubpd        %zmm3{cdab}, %zmm6, %zmm6{%k6}
	vmovnrngoapd  %zmm5,    0(%rbp)
	vmovnrngoapd  %zmm7,  128(%rbp)
	vmovnrngoapd  %zmm6,   64(%rbp)
	addq          $192, %rbx
	addq          $192, %rbp
	vmovapd          0(%rbx), %zmm2
	vmovapd        128(%rbx), %zmm4
	vmovapd         64(%rbx), %zmm3
	vmulpd        %zmm2, %zmm0, %zmm5
	vmulpd        %zmm2, %zmm1, %zmm2
	vaddpd        %zmm2{cdab}, %zmm5, %zmm5{%k7}
	vsubpd        %zmm2{cdab}, %zmm5, %zmm5{%k6}
	vmulpd        %zmm4, %zmm0, %zmm7
	vmulpd        %zmm4, %zmm1, %zmm4
	vaddpd        %zmm4{cdab}, %zmm7, %zmm7{%k7}
	vsubpd        %zmm4{cdab}, %zmm7, %zmm7{%k6}
	vmulpd        %zmm3, %zmm0, %zmm6
	vmulpd        %zmm3, %zmm1, %zmm3
	vaddpd        %zmm3{cdab}, %zmm6, %zmm6{%k7}
	vsubpd        %zmm3{cdab}, %zmm6, %zmm6{%k6}
	vmovnrngoapd  %zmm5,    0(%rbp)
	vmovnrngoapd  %zmm7,  128(%rbp)
	vmovnrngoapd  %zmm6,   64(%rbp)
	addq          $192, %rbx
	addq          $192, %rbp
	addq          $-1, %r10
	orq           %r10, %r10
	jg     _vmx_tmass_lab2
_vmx_tmass_lab1:
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
.type vmx_tmass,@function
.size vmx_tmass,.-vmx_tmass
.data
