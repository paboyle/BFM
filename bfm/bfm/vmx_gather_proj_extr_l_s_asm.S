/* 
 * BAGEL machine generated output.  
 *   
 * It is provided under the GNU pubic License V2  
 * It is provided as is, and is not guaranteed fit for any purpose.
 * BAGEL was written by Peter Boyle  
 */  
.text
.align 16,0x90
.globl vmx_gather_proj_extr_l_s
vmx_gather_proj_extr_l_s:
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
	addq          $0, %rsp
	movq          %rdi, %rcx
	movq             0(%rcx), %r8
	addq          $8, %rcx
	movq             0(%rcx), %rbx
	addq          $8, %rcx
	movq             0(%rcx), %r9
	addq          $8, %rcx
	movq             0(%rcx), %rbp
	addq          $8, %rcx
	movq             0(%rcx), %rax
	addq          $8, %rcx
	movq             0(%rcx), %r10
	addq          $8, %rcx
	movq             0(%rcx), %r13
	addq          $8, %rcx
	movq             0(%rcx), %r12
	addq          $8, %rcx
	movq             0(%rcx), %r14
	addq          $8, %rcx
	vmovapd          0(%r12), %zmm18
	orq           %rax, %rax
	jle    _vmx_gather_proj_extr_l_s_lab0
	imulq         $96, %r9, %r9
	imulq         %r10, %r9
	addq          %r9, %rbp
	addq          %r9, %r14
	movq          %r13, %r15
	orq           %r15, %r15
	je     _vmx_gather_proj_extr_l_s_lab2
	addq          $-1, %r15
	orq           %r15, %r15
	je     _vmx_gather_proj_extr_l_s_lab3
	addq          $-1, %r15
	orq           %r15, %r15
	je     _vmx_gather_proj_extr_l_s_lab4
	addq          $-1, %r15
	orq           %r15, %r15
	je     _vmx_gather_proj_extr_l_s_lab5
	addq          $-1, %r15
	orq           %r15, %r15
	je     _vmx_gather_proj_extr_l_s_lab6
	addq          $-1, %r15
	orq           %r15, %r15
	je     _vmx_gather_proj_extr_l_s_lab7
	addq          $-1, %r15
	orq           %r15, %r15
	je     _vmx_gather_proj_extr_l_s_lab8
	addq          $-1, %r15
	orq           %r15, %r15
	je     _vmx_gather_proj_extr_l_s_lab9
	addq          $-1, %r15
	jmp    _vmx_gather_proj_extr_l_s_lab1
_vmx_gather_proj_extr_l_s_lab2:
_vmx_gather_proj_extr_l_s_lab10:
	movq             0(%r8), %r11
	imulq         %r10, %r11
	imulq         $384, %r11, %r11
	addq          %rbx, %r11
	addq          $8, %r8
	movq          %r10, %r15
_vmx_gather_proj_extr_l_s_lab11:
	vmovaps          0(%r11), %zmm0
	vmovaps         64(%r11), %zmm2
	vmovaps        128(%r11), %zmm4
	vmovaps        192(%r11), %zmm6
	vmovaps        256(%r11), %zmm8
	vmovaps        320(%r11), %zmm10
	vmovaps         32(%r11), %zmm1
	vmovaps         96(%r11), %zmm3
	vmovaps        160(%r11), %zmm5
	vmovaps        224(%r11), %zmm7
	vmovaps        288(%r11), %zmm9
	vmovaps        352(%r11), %zmm11
	addq          $384, %r11
	vsubpd        %zmm6{cdab}, %zmm3, %zmm15{%k7}
	vaddpd        %zmm6{cdab}, %zmm3, %zmm15{%k6}
	vsubpd        %zmm7{cdab}, %zmm4, %zmm16{%k7}
	vaddpd        %zmm7{cdab}, %zmm4, %zmm16{%k6}
	vsubpd        %zmm8{cdab}, %zmm5, %zmm17{%k7}
	vaddpd        %zmm8{cdab}, %zmm5, %zmm17{%k6}
	vsubpd        %zmm9{cdab}, %zmm0, %zmm12{%k7}
	vaddpd        %zmm9{cdab}, %zmm0, %zmm12{%k6}
	vsubpd        %zmm10{cdab}, %zmm1, %zmm13{%k7}
	vaddpd        %zmm10{cdab}, %zmm1, %zmm13{%k6}
	vsubpd        %zmm11{cdab}, %zmm2, %zmm14{%k7}
	vaddpd        %zmm11{cdab}, %zmm2, %zmm14{%k6}
	vblendmpd     %zmm13{badc}, %zmm13, %zmm0
	vblendmpd     %zmm12, %zmm0, %zmm0{%k4}
	vblendmpd     %zmm15{badc}, %zmm15, %zmm1
	vblendmpd     %zmm14, %zmm1, %zmm1{%k4}
	vblendmpd     %zmm17{badc}, %zmm17, %zmm2
	vblendmpd     %zmm16, %zmm2, %zmm2{%k4}
	vmovnrngoaps  %zmm0,    0(%rbp)
	vmovnrngoaps  %zmm1,   32(%rbp)
	vmovnrngoaps  %zmm2,   64(%rbp)
	vblendmpd     %zmm12{badc}, %zmm12, %zmm0
	vblendmpd     %zmm13, %zmm0, %zmm0{%k4}
	vblendmpd     %zmm14{badc}, %zmm14, %zmm1
	vblendmpd     %zmm15, %zmm1, %zmm1{%k4}
	vblendmpd     %zmm16{badc}, %zmm16, %zmm2
	vblendmpd     %zmm17, %zmm2, %zmm2{%k4}
	vmovnrngoaps  %zmm0,    0(%r14)
	vmovnrngoaps  %zmm1,   32(%r14)
	vmovnrngoaps  %zmm2,   64(%r14)
	addq          $96, %rbp
	addq          $96, %r14
	addq          $-1, %r15
	orq           %r15, %r15
	jg     _vmx_gather_proj_extr_l_s_lab11
	addq          $-1, %rax
	orq           %rax, %rax
	jg     _vmx_gather_proj_extr_l_s_lab10
	jmp    _vmx_gather_proj_extr_l_s_lab1
_vmx_gather_proj_extr_l_s_lab3:
_vmx_gather_proj_extr_l_s_lab12:
	movq             0(%r8), %r11
	imulq         %r10, %r11
	imulq         $384, %r11, %r11
	addq          %rbx, %r11
	addq          $8, %r8
	movq          %r10, %r15
_vmx_gather_proj_extr_l_s_lab13:
	vmovaps          0(%r11), %zmm0
	vmovaps         64(%r11), %zmm2
	vmovaps        128(%r11), %zmm4
	vmovaps        192(%r11), %zmm6
	vmovaps        256(%r11), %zmm8
	vmovaps        320(%r11), %zmm10
	vmovaps         32(%r11), %zmm1
	vmovaps         96(%r11), %zmm3
	vmovaps        160(%r11), %zmm5
	vmovaps        224(%r11), %zmm7
	vmovaps        288(%r11), %zmm9
	vmovaps        352(%r11), %zmm11
	addq          $384, %r11
	vaddpd        %zmm6{cdab}, %zmm3, %zmm15{%k7}
	vsubpd        %zmm6{cdab}, %zmm3, %zmm15{%k6}
	vaddpd        %zmm7{cdab}, %zmm4, %zmm16{%k7}
	vsubpd        %zmm7{cdab}, %zmm4, %zmm16{%k6}
	vaddpd        %zmm8{cdab}, %zmm5, %zmm17{%k7}
	vsubpd        %zmm8{cdab}, %zmm5, %zmm17{%k6}
	vaddpd        %zmm9{cdab}, %zmm0, %zmm12{%k7}
	vsubpd        %zmm9{cdab}, %zmm0, %zmm12{%k6}
	vaddpd        %zmm10{cdab}, %zmm1, %zmm13{%k7}
	vsubpd        %zmm10{cdab}, %zmm1, %zmm13{%k6}
	vaddpd        %zmm11{cdab}, %zmm2, %zmm14{%k7}
	vsubpd        %zmm11{cdab}, %zmm2, %zmm14{%k6}
	vblendmpd     %zmm13{badc}, %zmm13, %zmm0
	vblendmpd     %zmm12, %zmm0, %zmm0{%k4}
	vblendmpd     %zmm15{badc}, %zmm15, %zmm1
	vblendmpd     %zmm14, %zmm1, %zmm1{%k4}
	vblendmpd     %zmm17{badc}, %zmm17, %zmm2
	vblendmpd     %zmm16, %zmm2, %zmm2{%k4}
	vmovnrngoaps  %zmm0,    0(%rbp)
	vmovnrngoaps  %zmm1,   32(%rbp)
	vmovnrngoaps  %zmm2,   64(%rbp)
	vblendmpd     %zmm12{badc}, %zmm12, %zmm0
	vblendmpd     %zmm13, %zmm0, %zmm0{%k4}
	vblendmpd     %zmm14{badc}, %zmm14, %zmm1
	vblendmpd     %zmm15, %zmm1, %zmm1{%k4}
	vblendmpd     %zmm16{badc}, %zmm16, %zmm2
	vblendmpd     %zmm17, %zmm2, %zmm2{%k4}
	vmovnrngoaps  %zmm0,    0(%r14)
	vmovnrngoaps  %zmm1,   32(%r14)
	vmovnrngoaps  %zmm2,   64(%r14)
	addq          $96, %rbp
	addq          $96, %r14
	addq          $-1, %r15
	orq           %r15, %r15
	jg     _vmx_gather_proj_extr_l_s_lab13
	addq          $-1, %rax
	orq           %rax, %rax
	jg     _vmx_gather_proj_extr_l_s_lab12
	jmp    _vmx_gather_proj_extr_l_s_lab1
_vmx_gather_proj_extr_l_s_lab4:
_vmx_gather_proj_extr_l_s_lab14:
	movq             0(%r8), %r11
	imulq         %r10, %r11
	imulq         $384, %r11, %r11
	addq          %rbx, %r11
	addq          $8, %r8
	movq          %r10, %r15
_vmx_gather_proj_extr_l_s_lab15:
	vmovaps          0(%r11), %zmm0
	vmovaps         64(%r11), %zmm2
	vmovaps        128(%r11), %zmm4
	vmovaps        192(%r11), %zmm6
	vmovaps        256(%r11), %zmm8
	vmovaps        320(%r11), %zmm10
	vmovaps         32(%r11), %zmm1
	vmovaps         96(%r11), %zmm3
	vmovaps        160(%r11), %zmm5
	vmovaps        224(%r11), %zmm7
	vmovaps        288(%r11), %zmm9
	vmovaps        352(%r11), %zmm11
	addq          $384, %r11
	vaddpd        %zmm6, %zmm3, %zmm15
	vaddpd        %zmm7, %zmm4, %zmm16
	vaddpd        %zmm8, %zmm5, %zmm17
	vsubpd        %zmm9, %zmm0, %zmm12
	vsubpd        %zmm10, %zmm1, %zmm13
	vsubpd        %zmm11, %zmm2, %zmm14
	vblendmpd     %zmm13{badc}, %zmm13, %zmm0
	vblendmpd     %zmm12, %zmm0, %zmm0{%k4}
	vblendmpd     %zmm15{badc}, %zmm15, %zmm1
	vblendmpd     %zmm14, %zmm1, %zmm1{%k4}
	vblendmpd     %zmm17{badc}, %zmm17, %zmm2
	vblendmpd     %zmm16, %zmm2, %zmm2{%k4}
	vmovnrngoaps  %zmm0,    0(%rbp)
	vmovnrngoaps  %zmm1,   32(%rbp)
	vmovnrngoaps  %zmm2,   64(%rbp)
	vblendmpd     %zmm12{badc}, %zmm12, %zmm0
	vblendmpd     %zmm13, %zmm0, %zmm0{%k4}
	vblendmpd     %zmm14{badc}, %zmm14, %zmm1
	vblendmpd     %zmm15, %zmm1, %zmm1{%k4}
	vblendmpd     %zmm16{badc}, %zmm16, %zmm2
	vblendmpd     %zmm17, %zmm2, %zmm2{%k4}
	vmovnrngoaps  %zmm0,    0(%r14)
	vmovnrngoaps  %zmm1,   32(%r14)
	vmovnrngoaps  %zmm2,   64(%r14)
	addq          $96, %rbp
	addq          $96, %r14
	addq          $-1, %r15
	orq           %r15, %r15
	jg     _vmx_gather_proj_extr_l_s_lab15
	addq          $-1, %rax
	orq           %rax, %rax
	jg     _vmx_gather_proj_extr_l_s_lab14
	jmp    _vmx_gather_proj_extr_l_s_lab1
_vmx_gather_proj_extr_l_s_lab5:
_vmx_gather_proj_extr_l_s_lab16:
	movq             0(%r8), %r11
	imulq         %r10, %r11
	imulq         $384, %r11, %r11
	addq          %rbx, %r11
	addq          $8, %r8
	movq          %r10, %r15
_vmx_gather_proj_extr_l_s_lab17:
	vmovaps          0(%r11), %zmm0
	vmovaps         64(%r11), %zmm2
	vmovaps        128(%r11), %zmm4
	vmovaps        192(%r11), %zmm6
	vmovaps        256(%r11), %zmm8
	vmovaps        320(%r11), %zmm10
	vmovaps         32(%r11), %zmm1
	vmovaps         96(%r11), %zmm3
	vmovaps        160(%r11), %zmm5
	vmovaps        224(%r11), %zmm7
	vmovaps        288(%r11), %zmm9
	vmovaps        352(%r11), %zmm11
	addq          $384, %r11
	vsubpd        %zmm6, %zmm3, %zmm15
	vsubpd        %zmm7, %zmm4, %zmm16
	vsubpd        %zmm8, %zmm5, %zmm17
	vaddpd        %zmm9, %zmm0, %zmm12
	vaddpd        %zmm10, %zmm1, %zmm13
	vaddpd        %zmm11, %zmm2, %zmm14
	vblendmpd     %zmm13{badc}, %zmm13, %zmm0
	vblendmpd     %zmm12, %zmm0, %zmm0{%k4}
	vblendmpd     %zmm15{badc}, %zmm15, %zmm1
	vblendmpd     %zmm14, %zmm1, %zmm1{%k4}
	vblendmpd     %zmm17{badc}, %zmm17, %zmm2
	vblendmpd     %zmm16, %zmm2, %zmm2{%k4}
	vmovnrngoaps  %zmm0,    0(%rbp)
	vmovnrngoaps  %zmm1,   32(%rbp)
	vmovnrngoaps  %zmm2,   64(%rbp)
	vblendmpd     %zmm12{badc}, %zmm12, %zmm0
	vblendmpd     %zmm13, %zmm0, %zmm0{%k4}
	vblendmpd     %zmm14{badc}, %zmm14, %zmm1
	vblendmpd     %zmm15, %zmm1, %zmm1{%k4}
	vblendmpd     %zmm16{badc}, %zmm16, %zmm2
	vblendmpd     %zmm17, %zmm2, %zmm2{%k4}
	vmovnrngoaps  %zmm0,    0(%r14)
	vmovnrngoaps  %zmm1,   32(%r14)
	vmovnrngoaps  %zmm2,   64(%r14)
	addq          $96, %rbp
	addq          $96, %r14
	addq          $-1, %r15
	orq           %r15, %r15
	jg     _vmx_gather_proj_extr_l_s_lab17
	addq          $-1, %rax
	orq           %rax, %rax
	jg     _vmx_gather_proj_extr_l_s_lab16
	jmp    _vmx_gather_proj_extr_l_s_lab1
_vmx_gather_proj_extr_l_s_lab6:
_vmx_gather_proj_extr_l_s_lab18:
	movq             0(%r8), %r11
	imulq         %r10, %r11
	imulq         $384, %r11, %r11
	addq          %rbx, %r11
	addq          $8, %r8
	movq          %r10, %r15
_vmx_gather_proj_extr_l_s_lab19:
	vmovaps          0(%r11), %zmm0
	vmovaps         64(%r11), %zmm2
	vmovaps        128(%r11), %zmm4
	vmovaps        192(%r11), %zmm6
	vmovaps        256(%r11), %zmm8
	vmovaps        320(%r11), %zmm10
	vmovaps         32(%r11), %zmm1
	vmovaps         96(%r11), %zmm3
	vmovaps        160(%r11), %zmm5
	vmovaps        224(%r11), %zmm7
	vmovaps        288(%r11), %zmm9
	vmovaps        352(%r11), %zmm11
	addq          $384, %r11
	vsubpd        %zmm6{cdab}, %zmm0, %zmm12{%k7}
	vaddpd        %zmm6{cdab}, %zmm0, %zmm12{%k6}
	vsubpd        %zmm7{cdab}, %zmm1, %zmm13{%k7}
	vaddpd        %zmm7{cdab}, %zmm1, %zmm13{%k6}
	vsubpd        %zmm8{cdab}, %zmm2, %zmm14{%k7}
	vaddpd        %zmm8{cdab}, %zmm2, %zmm14{%k6}
	vaddpd        %zmm9{cdab}, %zmm3, %zmm15{%k7}
	vsubpd        %zmm9{cdab}, %zmm3, %zmm15{%k6}
	vaddpd        %zmm10{cdab}, %zmm4, %zmm16{%k7}
	vsubpd        %zmm10{cdab}, %zmm4, %zmm16{%k6}
	vaddpd        %zmm11{cdab}, %zmm5, %zmm17{%k7}
	vsubpd        %zmm11{cdab}, %zmm5, %zmm17{%k6}
	vpermf32x4    $177, %zmm13, %zmm0
	vblendmpd     %zmm12, %zmm0, %zmm0{%k4}
	vpermf32x4    $177, %zmm15, %zmm1
	vblendmpd     %zmm14, %zmm1, %zmm1{%k4}
	vpermf32x4    $177, %zmm17, %zmm2
	vblendmpd     %zmm16, %zmm2, %zmm2{%k4}
	vmovnrngoaps  %zmm0,    0(%rbp)
	vmovnrngoaps  %zmm1,   32(%rbp)
	vmovnrngoaps  %zmm2,   64(%rbp)
	vpermf32x4    $177, %zmm12, %zmm0
	vblendmpd     %zmm13, %zmm0, %zmm0{%k4}
	vpermf32x4    $177, %zmm14, %zmm1
	vblendmpd     %zmm15, %zmm1, %zmm1{%k4}
	vpermf32x4    $177, %zmm16, %zmm2
	vblendmpd     %zmm17, %zmm2, %zmm2{%k4}
	vmovnrngoaps  %zmm0,    0(%r14)
	vmovnrngoaps  %zmm1,   32(%r14)
	vmovnrngoaps  %zmm2,   64(%r14)
	addq          $96, %rbp
	addq          $96, %r14
	addq          $-1, %r15
	orq           %r15, %r15
	jg     _vmx_gather_proj_extr_l_s_lab19
	addq          $-1, %rax
	orq           %rax, %rax
	jg     _vmx_gather_proj_extr_l_s_lab18
	jmp    _vmx_gather_proj_extr_l_s_lab1
_vmx_gather_proj_extr_l_s_lab7:
_vmx_gather_proj_extr_l_s_lab20:
	movq             0(%r8), %r11
	imulq         %r10, %r11
	imulq         $384, %r11, %r11
	addq          %rbx, %r11
	addq          $8, %r8
	movq          %r10, %r15
_vmx_gather_proj_extr_l_s_lab21:
	vmovaps          0(%r11), %zmm0
	vmovaps         64(%r11), %zmm2
	vmovaps        128(%r11), %zmm4
	vmovaps        192(%r11), %zmm6
	vmovaps        256(%r11), %zmm8
	vmovaps        320(%r11), %zmm10
	vmovaps         32(%r11), %zmm1
	vmovaps         96(%r11), %zmm3
	vmovaps        160(%r11), %zmm5
	vmovaps        224(%r11), %zmm7
	vmovaps        288(%r11), %zmm9
	vmovaps        352(%r11), %zmm11
	addq          $384, %r11
	vaddpd        %zmm6{cdab}, %zmm0, %zmm12{%k7}
	vsubpd        %zmm6{cdab}, %zmm0, %zmm12{%k6}
	vaddpd        %zmm7{cdab}, %zmm1, %zmm13{%k7}
	vsubpd        %zmm7{cdab}, %zmm1, %zmm13{%k6}
	vaddpd        %zmm8{cdab}, %zmm2, %zmm14{%k7}
	vsubpd        %zmm8{cdab}, %zmm2, %zmm14{%k6}
	vsubpd        %zmm9{cdab}, %zmm3, %zmm15{%k7}
	vaddpd        %zmm9{cdab}, %zmm3, %zmm15{%k6}
	vsubpd        %zmm10{cdab}, %zmm4, %zmm16{%k7}
	vaddpd        %zmm10{cdab}, %zmm4, %zmm16{%k6}
	vsubpd        %zmm11{cdab}, %zmm5, %zmm17{%k7}
	vaddpd        %zmm11{cdab}, %zmm5, %zmm17{%k6}
	vpermf32x4    $177, %zmm13, %zmm0
	vblendmpd     %zmm12, %zmm0, %zmm0{%k4}
	vpermf32x4    $177, %zmm15, %zmm1
	vblendmpd     %zmm14, %zmm1, %zmm1{%k4}
	vpermf32x4    $177, %zmm17, %zmm2
	vblendmpd     %zmm16, %zmm2, %zmm2{%k4}
	vmovnrngoaps  %zmm0,    0(%rbp)
	vmovnrngoaps  %zmm1,   32(%rbp)
	vmovnrngoaps  %zmm2,   64(%rbp)
	vpermf32x4    $177, %zmm12, %zmm0
	vblendmpd     %zmm13, %zmm0, %zmm0{%k4}
	vpermf32x4    $177, %zmm14, %zmm1
	vblendmpd     %zmm15, %zmm1, %zmm1{%k4}
	vpermf32x4    $177, %zmm16, %zmm2
	vblendmpd     %zmm17, %zmm2, %zmm2{%k4}
	vmovnrngoaps  %zmm0,    0(%r14)
	vmovnrngoaps  %zmm1,   32(%r14)
	vmovnrngoaps  %zmm2,   64(%r14)
	addq          $96, %rbp
	addq          $96, %r14
	addq          $-1, %r15
	orq           %r15, %r15
	jg     _vmx_gather_proj_extr_l_s_lab21
	addq          $-1, %rax
	orq           %rax, %rax
	jg     _vmx_gather_proj_extr_l_s_lab20
	jmp    _vmx_gather_proj_extr_l_s_lab1
_vmx_gather_proj_extr_l_s_lab8:
_vmx_gather_proj_extr_l_s_lab22:
	movq             0(%r8), %r11
	imulq         %r10, %r11
	imulq         $384, %r11, %r11
	addq          %rbx, %r11
	addq          $8, %r8
	movq          %r10, %r15
_vmx_gather_proj_extr_l_s_lab23:
	vmovaps          0(%r11), %zmm0
	vmovaps         64(%r11), %zmm2
	vmovaps        128(%r11), %zmm4
	vmovaps        192(%r11), %zmm6
	vmovaps        256(%r11), %zmm8
	vmovaps        320(%r11), %zmm10
	vmovaps         32(%r11), %zmm1
	vmovaps         96(%r11), %zmm3
	vmovaps        160(%r11), %zmm5
	vmovaps        224(%r11), %zmm7
	vmovaps        288(%r11), %zmm9
	vmovaps        352(%r11), %zmm11
	addq          $384, %r11
	vaddpd        %zmm6, %zmm0, %zmm12
	vaddpd        %zmm7, %zmm1, %zmm13
	vaddpd        %zmm8, %zmm2, %zmm14
	vaddpd        %zmm9, %zmm3, %zmm15
	vaddpd        %zmm10, %zmm4, %zmm16
	vaddpd        %zmm11, %zmm5, %zmm17
	vpermf32x4    $78, %zmm13, %zmm0
	vblendmpd     %zmm12, %zmm0, %zmm0{%k2}
	vpermf32x4    $78, %zmm15, %zmm1
	vblendmpd     %zmm14, %zmm1, %zmm1{%k2}
	vpermf32x4    $78, %zmm17, %zmm2
	vblendmpd     %zmm16, %zmm2, %zmm2{%k2}
	vmovnrngoaps  %zmm0,    0(%rbp)
	vmovnrngoaps  %zmm1,   32(%rbp)
	vmovnrngoaps  %zmm2,   64(%rbp)
	vpermf32x4    $78, %zmm12, %zmm0
	vblendmpd     %zmm13, %zmm0, %zmm0{%k2}
	vpermf32x4    $78, %zmm14, %zmm1
	vblendmpd     %zmm15, %zmm1, %zmm1{%k2}
	vpermf32x4    $78, %zmm16, %zmm2
	vblendmpd     %zmm17, %zmm2, %zmm2{%k2}
	vmovnrngoaps  %zmm0,    0(%r14)
	vmovnrngoaps  %zmm1,   32(%r14)
	vmovnrngoaps  %zmm2,   64(%r14)
	addq          $96, %rbp
	addq          $96, %r14
	addq          $-1, %r15
	orq           %r15, %r15
	jg     _vmx_gather_proj_extr_l_s_lab23
	addq          $-1, %rax
	orq           %rax, %rax
	jg     _vmx_gather_proj_extr_l_s_lab22
	jmp    _vmx_gather_proj_extr_l_s_lab1
_vmx_gather_proj_extr_l_s_lab9:
_vmx_gather_proj_extr_l_s_lab24:
	movq             0(%r8), %r11
	imulq         %r10, %r11
	imulq         $384, %r11, %r11
	addq          %rbx, %r11
	addq          $8, %r8
	movq          %r10, %r15
_vmx_gather_proj_extr_l_s_lab25:
	vmovaps          0(%r11), %zmm0
	vmovaps         64(%r11), %zmm2
	vmovaps        128(%r11), %zmm4
	vmovaps        192(%r11), %zmm6
	vmovaps        256(%r11), %zmm8
	vmovaps        320(%r11), %zmm10
	vmovaps         32(%r11), %zmm1
	vmovaps         96(%r11), %zmm3
	vmovaps        160(%r11), %zmm5
	vmovaps        224(%r11), %zmm7
	vmovaps        288(%r11), %zmm9
	vmovaps        352(%r11), %zmm11
	addq          $384, %r11
	vsubpd        %zmm6, %zmm0, %zmm12
	vsubpd        %zmm7, %zmm1, %zmm13
	vsubpd        %zmm8, %zmm2, %zmm14
	vsubpd        %zmm9, %zmm3, %zmm15
	vsubpd        %zmm10, %zmm4, %zmm16
	vsubpd        %zmm11, %zmm5, %zmm17
	vpermf32x4    $78, %zmm13, %zmm0
	vblendmpd     %zmm12, %zmm0, %zmm0{%k2}
	vpermf32x4    $78, %zmm15, %zmm1
	vblendmpd     %zmm14, %zmm1, %zmm1{%k2}
	vpermf32x4    $78, %zmm17, %zmm2
	vblendmpd     %zmm16, %zmm2, %zmm2{%k2}
	vmovnrngoaps  %zmm0,    0(%rbp)
	vmovnrngoaps  %zmm1,   32(%rbp)
	vmovnrngoaps  %zmm2,   64(%rbp)
	vpermf32x4    $78, %zmm12, %zmm0
	vblendmpd     %zmm13, %zmm0, %zmm0{%k2}
	vpermf32x4    $78, %zmm14, %zmm1
	vblendmpd     %zmm15, %zmm1, %zmm1{%k2}
	vpermf32x4    $78, %zmm16, %zmm2
	vblendmpd     %zmm17, %zmm2, %zmm2{%k2}
	vmovnrngoaps  %zmm0,    0(%r14)
	vmovnrngoaps  %zmm1,   32(%r14)
	vmovnrngoaps  %zmm2,   64(%r14)
	addq          $96, %rbp
	addq          $96, %r14
	addq          $-1, %r15
	orq           %r15, %r15
	jg     _vmx_gather_proj_extr_l_s_lab25
	addq          $-1, %rax
	orq           %rax, %rax
	jg     _vmx_gather_proj_extr_l_s_lab24
	jmp    _vmx_gather_proj_extr_l_s_lab1
_vmx_gather_proj_extr_l_s_lab1:
_vmx_gather_proj_extr_l_s_lab0:
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
.type vmx_gather_proj_extr_l_s,@function
.size vmx_gather_proj_extr_l_s,.-vmx_gather_proj_extr_l_s
.data
