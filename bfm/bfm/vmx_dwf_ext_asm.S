/*Exterior links*/
/*Add to result field*/
/* 
 * BAGEL machine generated output.  
 *   
 * It is provided under the GNU pubic License V2  
 * It is provided as is, and is not guaranteed fit for any purpose.
 * BAGEL was written by Peter Boyle  
 */  
.text
.align 16,0x90
.globl vmx_dwf_ext
vmx_dwf_ext:
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

	addq          $-2688, %rsp
	movq          %rbx, 2624(%rsp)
	movq          %rbp, 2632(%rsp)
	movq          %r12, 2656(%rsp)
	movq          %r13, 2664(%rsp)
	movq          %r14, 2672(%rsp)
	movq          %r15, 2680(%rsp)
	addq          $0, %rsp
	movq          %rdi, %rbp
	addq          $8, %rbp
	movq             0(%rbp), %rdx
	addq          $8, %rbp
	movq             0(%rbp), %rsi
	addq          $8, %rbp
	movq          %rsi,   40(%rsp)
	movq             0(%rbp), %r13
	addq          $8, %rbp
	movq             0(%rbp), %r9
	addq          $8, %rbp
	movq             0(%rbp), %rcx
	addq          $8, %rbp
	movq          %rcx,   32(%rsp)
	movq             0(%rbp), %r8
	addq          $8, %rbp
	vmovapd          0(%r8), %zmm27
	movq             0(%rbp), %rdi
	addq          $8, %rbp
	movq          %rdi,   96(%rsp)
	movq             0(%rbp), %rbx
	addq          $8, %rbp
	vmovapd        128(%r8), %zmm24
	movq             0(%rbp), %r15
	orq           %r13, %r13
	movq          %r15,   48(%rsp)
	movq          %r15,   80(%rsp)
	jle    _vmx_dwf_ext_lab0
_vmx_dwf_ext_lab1:
	movq          %r13,   64(%rsp)
	movq            96(%rsp), %r11
	movq            40(%rsp), %rsi
	movq             0(%r11), %rdi
	addq          $8, %r11
	movq            32(%rsp), %rcx
	movq          %r11,   96(%rsp)
	imulq         $576, %rdi, %r11
	movq            48(%rsp), %r15
	imulq         $8, %r11, %r11
	addq          %r11, %rsi
	movq          %rsi, %r11
	vprefetch1       0(%r11)
	vprefetch1      64(%r11)
	vprefetch1     128(%r11)
	vprefetch1     192(%r11)
	vprefetch1     256(%r11)
	vprefetch1     320(%r11)
	vprefetch1     384(%r11)
	vprefetch1     448(%r11)
	vprefetch1     512(%r11)
	addq          $576, %r11
	vprefetch1       0(%r11)
	vprefetch1      64(%r11)
	vprefetch1     128(%r11)
	vprefetch1     192(%r11)
	vprefetch1     256(%r11)
	vprefetch1     320(%r11)
	vprefetch1     384(%r11)
	vprefetch1     448(%r11)
	vprefetch1     512(%r11)
	addq          $576, %r11
	vprefetch1       0(%r11)
	vprefetch1      64(%r11)
	vprefetch1     128(%r11)
	vprefetch1     192(%r11)
	vprefetch1     256(%r11)
	vprefetch1     320(%r11)
	vprefetch1     384(%r11)
	vprefetch1     448(%r11)
	vprefetch1     512(%r11)
	addq          $576, %r11
	vprefetch1       0(%r11)
	vprefetch1      64(%r11)
	vprefetch1     128(%r11)
	vprefetch1     192(%r11)
	vprefetch1     256(%r11)
	vprefetch1     320(%r11)
	vprefetch1     384(%r11)
	vprefetch1     448(%r11)
	vprefetch1     512(%r11)
	addq          $576, %r11
	vprefetch1       0(%r11)
	vprefetch1      64(%r11)
	vprefetch1     128(%r11)
	vprefetch1     192(%r11)
	vprefetch1     256(%r11)
	vprefetch1     320(%r11)
	vprefetch1     384(%r11)
	vprefetch1     448(%r11)
	vprefetch1     512(%r11)
	addq          $576, %r11
	vprefetch1       0(%r11)
	vprefetch1      64(%r11)
	vprefetch1     128(%r11)
	vprefetch1     192(%r11)
	vprefetch1     256(%r11)
	vprefetch1     320(%r11)
	vprefetch1     384(%r11)
	vprefetch1     448(%r11)
	vprefetch1     512(%r11)
	addq          $576, %r11
	vprefetch1       0(%r11)
	vprefetch1      64(%r11)
	vprefetch1     128(%r11)
	vprefetch1     192(%r11)
	vprefetch1     256(%r11)
	vprefetch1     320(%r11)
	vprefetch1     384(%r11)
	vprefetch1     448(%r11)
	vprefetch1     512(%r11)
	addq          $576, %r11
	vprefetch1       0(%r11)
	vprefetch1      64(%r11)
	vprefetch1     128(%r11)
	vprefetch1     192(%r11)
	vprefetch1     256(%r11)
	vprefetch1     320(%r11)
	vprefetch1     384(%r11)
	vprefetch1     448(%r11)
	vprefetch1     512(%r11)
	addq          $576, %r11
	imulq         $128, %rdi, %r11
	addq          %r11, %rcx
	imulq         $768, %rdi, %r11
	imulq         %r9, %r11
	movq            80(%rcx), %rbp
	addq          %r11, %r15
	movq          $0,%r11
	movq          %r15,   80(%rsp)
	movq          %rbp, %r13
	andq          $1, %r13
	movq             0(%rcx), %r12
	movq          %rdx, %r15
	orq           %r13, %r13
	jle    _vmx_dwf_ext_lab2
	movq          %rbx, %r15
_vmx_dwf_ext_lab2:
	addq          %r12, %r15
	movq          $-1,%r14
	movq          %r15,  128(%rsp)
	orq           %r13, %r13
	jle    _vmx_dwf_ext_lab3
	movq          %r11, %r14
_vmx_dwf_ext_lab3:
	movq          %r14,  256(%rsp)
	movq          $384,%r14
	orq           %r13, %r13
	movq          %r14, %r15
	addq          %r14, %r15
	jle    _vmx_dwf_ext_lab4
	movq          %r14, %r15
_vmx_dwf_ext_lab4:
	movq          %r15,  192(%rsp)
	movq          %rbp, %r13
	andq          $2, %r13
	movq             8(%rcx), %r12
	movq          %rdx, %r15
	orq           %r13, %r13
	jle    _vmx_dwf_ext_lab5
	movq          %rbx, %r15
_vmx_dwf_ext_lab5:
	addq          %r12, %r15
	movq          $-1,%r14
	movq          %r15,  136(%rsp)
	orq           %r13, %r13
	jle    _vmx_dwf_ext_lab6
	movq          %r11, %r14
_vmx_dwf_ext_lab6:
	movq          %r14,  264(%rsp)
	movq          $384,%r14
	orq           %r13, %r13
	movq          %r14, %r15
	addq          %r14, %r15
	jle    _vmx_dwf_ext_lab7
	movq          %r14, %r15
_vmx_dwf_ext_lab7:
	movq          %r15,  200(%rsp)
	movq          %rbp, %r13
	andq          $4, %r13
	movq            16(%rcx), %r12
	movq          %rdx, %r15
	orq           %r13, %r13
	jle    _vmx_dwf_ext_lab8
	movq          %rbx, %r15
_vmx_dwf_ext_lab8:
	addq          %r12, %r15
	movq          $-1,%r14
	movq          %r15,  144(%rsp)
	orq           %r13, %r13
	jle    _vmx_dwf_ext_lab9
	movq          %r11, %r14
_vmx_dwf_ext_lab9:
	movq          %r14,  272(%rsp)
	movq          $384,%r14
	orq           %r13, %r13
	movq          %r14, %r15
	addq          %r14, %r15
	jle    _vmx_dwf_ext_lab10
	movq          %r14, %r15
_vmx_dwf_ext_lab10:
	movq          %r15,  208(%rsp)
	movq          %rbp, %r13
	andq          $8, %r13
	movq            24(%rcx), %r12
	movq          %rdx, %r15
	orq           %r13, %r13
	jle    _vmx_dwf_ext_lab11
	movq          %rbx, %r15
_vmx_dwf_ext_lab11:
	addq          %r12, %r15
	movq          $-1,%r14
	movq          %r15,  152(%rsp)
	orq           %r13, %r13
	jle    _vmx_dwf_ext_lab12
	movq          %r11, %r14
_vmx_dwf_ext_lab12:
	movq          %r14,  280(%rsp)
	movq          $384,%r14
	orq           %r13, %r13
	movq          %r14, %r15
	addq          %r14, %r15
	jle    _vmx_dwf_ext_lab13
	movq          %r14, %r15
_vmx_dwf_ext_lab13:
	movq          %r15,  216(%rsp)
	movq          %rbp, %r13
	andq          $16, %r13
	movq            32(%rcx), %r12
	movq          %rdx, %r15
	orq           %r13, %r13
	jle    _vmx_dwf_ext_lab14
	movq          %rbx, %r15
_vmx_dwf_ext_lab14:
	addq          %r12, %r15
	movq          $-1,%r14
	movq          %r15,  160(%rsp)
	orq           %r13, %r13
	jle    _vmx_dwf_ext_lab15
	movq          %r11, %r14
_vmx_dwf_ext_lab15:
	movq          %r14,  288(%rsp)
	movq          $384,%r14
	orq           %r13, %r13
	movq          %r14, %r15
	addq          %r14, %r15
	jle    _vmx_dwf_ext_lab16
	movq          %r14, %r15
_vmx_dwf_ext_lab16:
	movq          %r15,  224(%rsp)
	movq          %rbp, %r13
	andq          $32, %r13
	movq            40(%rcx), %r12
	movq          %rdx, %r15
	orq           %r13, %r13
	jle    _vmx_dwf_ext_lab17
	movq          %rbx, %r15
_vmx_dwf_ext_lab17:
	addq          %r12, %r15
	movq          $-1,%r14
	movq          %r15,  168(%rsp)
	orq           %r13, %r13
	jle    _vmx_dwf_ext_lab18
	movq          %r11, %r14
_vmx_dwf_ext_lab18:
	movq          %r14,  296(%rsp)
	movq          $384,%r14
	orq           %r13, %r13
	movq          %r14, %r15
	addq          %r14, %r15
	jle    _vmx_dwf_ext_lab19
	movq          %r14, %r15
_vmx_dwf_ext_lab19:
	movq          %r15,  232(%rsp)
	movq          %rbp, %r13
	andq          $64, %r13
	movq            48(%rcx), %r12
	movq          %rdx, %r15
	orq           %r13, %r13
	jle    _vmx_dwf_ext_lab20
	movq          %rbx, %r15
_vmx_dwf_ext_lab20:
	addq          %r12, %r15
	movq          $-1,%r14
	movq          %r15,  176(%rsp)
	orq           %r13, %r13
	jle    _vmx_dwf_ext_lab21
	movq          %r11, %r14
_vmx_dwf_ext_lab21:
	movq          %r14,  304(%rsp)
	movq          $384,%r14
	orq           %r13, %r13
	movq          %r14, %r15
	addq          %r14, %r15
	jle    _vmx_dwf_ext_lab22
	movq          %r14, %r15
_vmx_dwf_ext_lab22:
	movq          %r15,  240(%rsp)
	movq          %rbp, %r13
	andq          $128, %r13
	movq            56(%rcx), %r12
	movq          %rdx, %r15
	orq           %r13, %r13
	jle    _vmx_dwf_ext_lab23
	movq          %rbx, %r15
_vmx_dwf_ext_lab23:
	addq          %r12, %r15
	movq          $-1,%r14
	movq          %r15,  184(%rsp)
	orq           %r13, %r13
	jle    _vmx_dwf_ext_lab24
	movq          %r11, %r14
_vmx_dwf_ext_lab24:
	movq          %r14,  312(%rsp)
	movq          $384,%r14
	orq           %r13, %r13
	movq          %r14, %r15
	addq          %r14, %r15
	jle    _vmx_dwf_ext_lab25
	movq          %r14, %r15
_vmx_dwf_ext_lab25:
	movq          %r15,  248(%rsp)
	movq          %r9, %rax
_vmx_dwf_ext_lab26:
	vpxorq        %zmm6, %zmm6, %zmm6
	movq           128(%rsp), %r12
	vpxorq        %zmm7, %zmm7, %zmm7
	movq           192(%rsp), %r15
	vpxorq        %zmm8, %zmm8, %zmm8
	movq          %r9, %r14
	vpxorq        %zmm9, %zmm9, %zmm9
	subq          %rax, %r14
	vpxorq        %zmm10, %zmm10, %zmm10
	movq          %rsi, %r11
	vpxorq        %zmm11, %zmm11, %zmm11
	imulq         %r14, %r15
	vmovapd          0(%r11), %zmm24
	addq          %r15, %r12
	vmovapd         64(%r11), %zmm25
	orq           %r12, %r12
	vmovapd        128(%r11), %zmm26
	movq          %r12, %r10
	movq           136(%rsp), %r12
	vmovapd          0(%r10), %zmm12
	vmovapd         64(%r10), %zmm13
	vmovapd        128(%r10), %zmm14
	vmovapd        192(%r10), %zmm15
	vmovapd        256(%r10), %zmm16
	vmovapd        320(%r10), %zmm17
	vmovapd        384(%r10), %zmm18
	vmovapd        448(%r10), %zmm19
	vmovapd        512(%r10), %zmm20
	vpandd         256(%rsp){1to16} , %zmm18, %zmm18
	vmovapd        576(%r10), %zmm21
	vpandd         256(%rsp){1to16} , %zmm19, %zmm19
	vmovapd        640(%r10), %zmm22
	vpandd         256(%rsp){1to16} , %zmm20, %zmm20
	vmovapd        704(%r10), %zmm23
	vpandd         256(%rsp){1to16} , %zmm21, %zmm21
	movq           200(%rsp), %r15
	vpandd         256(%rsp){1to16} , %zmm22, %zmm22
	movq          %r9, %r14
	vpandd         256(%rsp){1to16} , %zmm23, %zmm23
	subq          %rax, %r14
	vsubpd        %zmm18{cdab}, %zmm15, %zmm15{%k7}
	movq          %r8, %r13
	vsubpd        %zmm21{cdab}, %zmm12, %zmm12{%k7}
	imulq         %r14, %r15
	vsubpd        %zmm22{cdab}, %zmm13, %zmm13{%k7}
	movq          %r9, %r14
	vsubpd        %zmm23{cdab}, %zmm14, %zmm14{%k7}
	addq          %r15, %r12
	vsubpd        %zmm19{cdab}, %zmm16, %zmm16{%k7}
	movq           208(%rsp), %r15
	vsubpd        %zmm20{cdab}, %zmm17, %zmm17{%k7}
	vprefetch0       0(%r12)
	vaddpd        %zmm18{cdab}, %zmm15, %zmm15{%k6}
	vprefetch1     768(%r12)
	vaddpd        %zmm21{cdab}, %zmm12, %zmm12{%k6}
	vprefetch0      64(%r12)
	vaddpd        %zmm22{cdab}, %zmm13, %zmm13{%k6}
	vprefetch1     832(%r12)
	vaddpd        %zmm23{cdab}, %zmm14, %zmm14{%k6}
	vprefetch0     128(%r12)
	vaddpd        %zmm19{cdab}, %zmm16, %zmm16{%k6}
	vprefetch1     896(%r12)
	vaddpd        %zmm20{cdab}, %zmm17, %zmm17{%k6}
	vprefetch0     192(%r12)
	vmulpd        %zmm15, %zmm24, %zmm21
	vprefetch1     960(%r12)
	vmulpd        %zmm12, %zmm24, %zmm18
	vprefetch0     256(%r12)
	vmulpd        %zmm12{cdab}, %zmm24, %zmm28
	vprefetch1    1024(%r12)
	vmulpd        %zmm15{cdab}, %zmm24, %zmm29
	vprefetch0     320(%r12)
	vfmadd231pd   %zmm13, %zmm25, %zmm18
	vmovapd        192(%r11), %zmm24
	vfmadd231pd        %zmm13{cdab}, %zmm25, %zmm28
	vprefetch1    1088(%r12)
	vfmadd231pd   %zmm16, %zmm25, %zmm21
	vprefetch0     384(%r12)
	vfmadd231pd        %zmm16{cdab}, %zmm25, %zmm29
	vprefetch1    1152(%r12)
	vfmadd231pd   %zmm14, %zmm26, %zmm18
	vmovapd        256(%r11), %zmm25
	vfmadd231pd        %zmm14{cdab}, %zmm26, %zmm28
	vprefetch0     448(%r12)
	vfmadd231pd   %zmm17, %zmm26, %zmm21
	vprefetch1    1216(%r12)
	vfmadd231pd        %zmm17{cdab}, %zmm26, %zmm29
	vprefetch0     512(%r12)
	vsubpd        %zmm18{cdab}, %zmm18, %zmm18{%k7}
	vmovapd        320(%r11), %zmm26
	vmulpd        %zmm12, %zmm24, %zmm19
	vprefetch1    1280(%r12)
	vsubpd        %zmm21{cdab}, %zmm21, %zmm21{%k7}
	vprefetch0     576(%r12)
	vaddpd        %zmm28{cdab}, %zmm28, %zmm18{%k6}
	vprefetch1    1344(%r12)
	vmulpd        %zmm12{cdab}, %zmm24, %zmm30
	vprefetch0     640(%r12)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm21{%k6}
	vprefetch1    1408(%r12)
	vmulpd        %zmm15, %zmm24, %zmm22
	vprefetch0     704(%r12)
	vmulpd        %zmm15{cdab}, %zmm24, %zmm31
	vprefetch1    1472(%r12)
	vfmadd231pd   %zmm13, %zmm25, %zmm19
	vmovapd        384(%r11), %zmm24
	vfmadd231pd        %zmm13{cdab}, %zmm25, %zmm30
	orq           %r12, %r12
	vfmadd231pd   %zmm16, %zmm25, %zmm22
	subq          %rax, %r14
	vfmadd231pd        %zmm16{cdab}, %zmm25, %zmm31
	movq          %r12, %r10
	vfmadd231pd   %zmm14, %zmm26, %zmm19
	vmovapd        448(%r11), %zmm25
	vfmadd231pd        %zmm14{cdab}, %zmm26, %zmm30
	movq           144(%rsp), %r12
	vfmadd231pd   %zmm17, %zmm26, %zmm22
	imulq         %r14, %r15
	vfmadd231pd        %zmm17{cdab}, %zmm26, %zmm31
	movq          %r9, %r14
	vsubpd        %zmm19{cdab}, %zmm19, %zmm19{%k7}
	vmovapd        512(%r11), %zmm26
	vmulpd        %zmm12, %zmm24, %zmm20
	addq          $576, %r11
	vsubpd        %zmm22{cdab}, %zmm22, %zmm22{%k7}
	addq          %r15, %r12
	vaddpd        %zmm30{cdab}, %zmm30, %zmm19{%k6}
	movq           216(%rsp), %r15
	vmulpd        %zmm12{cdab}, %zmm24, %zmm28
	vprefetch0       0(%r12)
	vaddpd        %zmm31{cdab}, %zmm31, %zmm22{%k6}
	vmovapd          0(%r10), %zmm12
	vmulpd        %zmm15, %zmm24, %zmm23
	vprefetch1     768(%r12)
	vmulpd        %zmm15{cdab}, %zmm24, %zmm29
	vprefetch0      64(%r12)
	vfmadd231pd   %zmm13, %zmm25, %zmm20
	vmovapd        192(%r10), %zmm15
	vfmadd231pd        %zmm13{cdab}, %zmm25, %zmm28
	vmovapd          0(%r11), %zmm24
	vfmadd231pd   %zmm16, %zmm25, %zmm23
	vmovapd         64(%r10), %zmm13
	vfmadd231pd        %zmm16{cdab}, %zmm25, %zmm29
	vprefetch1     832(%r12)
	vfmadd231pd   %zmm14, %zmm26, %zmm20
	vmovapd        256(%r10), %zmm16
	vfmadd231pd        %zmm14{cdab}, %zmm26, %zmm28
	vmovapd         64(%r11), %zmm25
	vfmadd231pd   %zmm17, %zmm26, %zmm23
	vmovapd        128(%r10), %zmm14
	vfmadd231pd        %zmm17{cdab}, %zmm26, %zmm29
	vprefetch0     128(%r12)
	vsubpd        %zmm20{cdab}, %zmm20, %zmm20{%k7}
	vmovapd        320(%r10), %zmm17
	vmovapd       %zmm18, %zmm0
	vmovapd        128(%r11), %zmm26
	vsubpd        %zmm23{cdab}, %zmm23, %zmm23{%k7}
	vprefetch1     896(%r12)
	vaddpd        %zmm28{cdab}, %zmm28, %zmm20{%k6}
	vprefetch0     192(%r12)
	vmovapd       %zmm19, %zmm1
	vprefetch1     960(%r12)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm23{%k6}
	vprefetch0     256(%r12)
	vmovapd       %zmm21, %zmm3
	vprefetch1    1024(%r12)
	vmovapd       %zmm22, %zmm4
	vprefetch0     320(%r12)
	vaddpd        %zmm21{cdab}, %zmm6, %zmm6{%k7}
	vprefetch1    1088(%r12)
	vmovapd       %zmm20, %zmm2
	vprefetch0     384(%r12)
	vaddpd        %zmm22{cdab}, %zmm7, %zmm7{%k7}
	vprefetch1    1152(%r12)
	vmovapd       %zmm23, %zmm5
	vprefetch0     448(%r12)
	vaddpd        %zmm23{cdab}, %zmm8, %zmm8{%k7}
	vprefetch1    1216(%r12)
	vaddpd        %zmm18{cdab}, %zmm9, %zmm9{%k7}
	vprefetch0     512(%r12)
	vsubpd        %zmm21{cdab}, %zmm6, %zmm6{%k6}
	vprefetch1    1280(%r12)
	vaddpd        %zmm19{cdab}, %zmm10, %zmm10{%k7}
	vmovapd        576(%r10), %zmm21
	vsubpd        %zmm22{cdab}, %zmm7, %zmm7{%k6}
	vprefetch0     576(%r12)
	vaddpd        %zmm20{cdab}, %zmm11, %zmm11{%k7}
	vmovapd        640(%r10), %zmm22
	vsubpd        %zmm23{cdab}, %zmm8, %zmm8{%k6}
	vprefetch1    1344(%r12)
	vsubpd        %zmm18{cdab}, %zmm9, %zmm9{%k6}
	vmovapd        704(%r10), %zmm23
	vpandd         264(%rsp){1to16} , %zmm21, %zmm21
	vmovapd        384(%r10), %zmm18
	vsubpd        %zmm19{cdab}, %zmm10, %zmm10{%k6}
	vprefetch0     640(%r12)
	vpandd         264(%rsp){1to16} , %zmm22, %zmm22
	vmovapd        448(%r10), %zmm19
	vsubpd        %zmm20{cdab}, %zmm11, %zmm11{%k6}
	vprefetch1    1408(%r12)
	vpandd         264(%rsp){1to16} , %zmm18, %zmm18
	vmovapd        512(%r10), %zmm20
	vpandd         264(%rsp){1to16} , %zmm19, %zmm19
	vprefetch0     704(%r12)
	vpandd         264(%rsp){1to16} , %zmm23, %zmm23
	vprefetch1    1472(%r12)
	vpandd         264(%rsp){1to16} , %zmm20, %zmm20
	orq           %r12, %r12
	vaddpd        %zmm21{cdab}, %zmm12, %zmm12{%k7}
	subq          %rax, %r14
	vaddpd        %zmm22{cdab}, %zmm13, %zmm13{%k7}
	movq          %r12, %r10
	vaddpd        %zmm23{cdab}, %zmm14, %zmm14{%k7}
	movq           152(%rsp), %r12
	vaddpd        %zmm18{cdab}, %zmm15, %zmm15{%k7}
	imulq         %r14, %r15
	vaddpd        %zmm19{cdab}, %zmm16, %zmm16{%k7}
	movq          %r9, %r14
	vaddpd        %zmm20{cdab}, %zmm17, %zmm17{%k7}
	addq          %r15, %r12
	vsubpd        %zmm21{cdab}, %zmm12, %zmm12{%k6}
	movq           224(%rsp), %r15
	vsubpd        %zmm22{cdab}, %zmm13, %zmm13{%k6}
	vprefetch0       0(%r12)
	vsubpd        %zmm23{cdab}, %zmm14, %zmm14{%k6}
	vprefetch1     768(%r12)
	vsubpd        %zmm18{cdab}, %zmm15, %zmm15{%k6}
	vprefetch0      64(%r12)
	vsubpd        %zmm19{cdab}, %zmm16, %zmm16{%k6}
	vprefetch1     832(%r12)
	vsubpd        %zmm20{cdab}, %zmm17, %zmm17{%k6}
	vprefetch0     128(%r12)
	vmulpd        %zmm12, %zmm24, %zmm18
	vprefetch1     896(%r12)
	vmulpd        %zmm12{cdab}, %zmm24, %zmm28
	vprefetch0     192(%r12)
	vfmadd231pd   %zmm13, %zmm25, %zmm18
	vprefetch1     960(%r12)
	vmulpd        %zmm15, %zmm24, %zmm21
	vprefetch0     256(%r12)
	vmulpd        %zmm15{cdab}, %zmm24, %zmm29
	vprefetch1    1024(%r12)
	vfmadd231pd        %zmm13{cdab}, %zmm25, %zmm28
	vmovapd        192(%r11), %zmm24
	vfmadd231pd   %zmm16, %zmm25, %zmm21
	vprefetch0     320(%r12)
	vfmadd231pd        %zmm16{cdab}, %zmm25, %zmm29
	vprefetch1    1088(%r12)
	vfmadd231pd   %zmm14, %zmm26, %zmm18
	vmovapd        256(%r11), %zmm25
	vfmadd231pd        %zmm14{cdab}, %zmm26, %zmm28
	vprefetch0     384(%r12)
	vfmadd231pd   %zmm17, %zmm26, %zmm21
	vprefetch1    1152(%r12)
	vfmadd231pd        %zmm17{cdab}, %zmm26, %zmm29
	vprefetch0     448(%r12)
	vsubpd        %zmm18{cdab}, %zmm18, %zmm18{%k7}
	vmovapd        320(%r11), %zmm26
	vmulpd        %zmm12, %zmm24, %zmm19
	vprefetch1    1216(%r12)
	vsubpd        %zmm21{cdab}, %zmm21, %zmm21{%k7}
	vprefetch0     512(%r12)
	vaddpd        %zmm28{cdab}, %zmm28, %zmm18{%k6}
	vprefetch1    1280(%r12)
	vmulpd        %zmm12{cdab}, %zmm24, %zmm28
	vprefetch0     576(%r12)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm21{%k6}
	vprefetch1    1344(%r12)
	vmulpd        %zmm15, %zmm24, %zmm22
	vprefetch0     640(%r12)
	vmulpd        %zmm15{cdab}, %zmm24, %zmm29
	vprefetch1    1408(%r12)
	vfmadd231pd   %zmm13, %zmm25, %zmm19
	vmovapd        384(%r11), %zmm24
	vfmadd231pd        %zmm13{cdab}, %zmm25, %zmm28
	vprefetch0     704(%r12)
	vfmadd231pd   %zmm16, %zmm25, %zmm22
	vprefetch1    1472(%r12)
	vfmadd231pd        %zmm16{cdab}, %zmm25, %zmm29
	orq           %r12, %r12
	vfmadd231pd   %zmm14, %zmm26, %zmm19
	vmovapd        448(%r11), %zmm25
	vfmadd231pd        %zmm14{cdab}, %zmm26, %zmm28
	subq          %rax, %r14
	vfmadd231pd   %zmm17, %zmm26, %zmm22
	addq          $448, %r13
	vfmadd231pd        %zmm17{cdab}, %zmm26, %zmm29
	imulq         %r14, %r15
	vsubpd        %zmm19{cdab}, %zmm19, %zmm19{%k7}
	vmovapd        512(%r11), %zmm26
	vmulpd        %zmm12, %zmm24, %zmm20
	addq          $576, %r11
	vsubpd        %zmm22{cdab}, %zmm22, %zmm22{%k7}
	movq          %r9, %r14
	vaddpd        %zmm28{cdab}, %zmm28, %zmm19{%k6}
	subq          %rax, %r14
	vmulpd        %zmm12{cdab}, %zmm24, %zmm28
	vaddpd        %zmm29{cdab}, %zmm29, %zmm22{%k6}
	vmovapd          0(%r10), %zmm12
	vmulpd        %zmm15, %zmm24, %zmm23
	vmulpd        %zmm15{cdab}, %zmm24, %zmm29
	vfmadd231pd   %zmm13, %zmm25, %zmm20
	vmovapd        192(%r10), %zmm15
	vfmadd231pd        %zmm13{cdab}, %zmm25, %zmm28
	vmovapd          0(%r11), %zmm24
	vfmadd231pd   %zmm16, %zmm25, %zmm23
	vmovapd         64(%r10), %zmm13
	vfmadd231pd        %zmm16{cdab}, %zmm25, %zmm29
	vfmadd231pd   %zmm14, %zmm26, %zmm20
	vmovapd        256(%r10), %zmm16
	vfmadd231pd        %zmm14{cdab}, %zmm26, %zmm28
	vmovapd         64(%r11), %zmm25
	vfmadd231pd   %zmm17, %zmm26, %zmm23
	vmovapd        128(%r10), %zmm14
	vfmadd231pd        %zmm17{cdab}, %zmm26, %zmm29
	vsubpd        %zmm20{cdab}, %zmm20, %zmm20{%k7}
	vmovapd        320(%r10), %zmm17
	vaddpd        %zmm18, %zmm0, %zmm0
	vmovapd        128(%r11), %zmm26
	vsubpd        %zmm23{cdab}, %zmm23, %zmm23{%k7}
	vaddpd        %zmm28{cdab}, %zmm28, %zmm20{%k6}
	vaddpd        %zmm19, %zmm1, %zmm1
	vaddpd        %zmm29{cdab}, %zmm29, %zmm23{%k6}
	vaddpd        %zmm21, %zmm3, %zmm3
	vaddpd        %zmm22, %zmm4, %zmm4
	vsubpd        %zmm21{cdab}, %zmm6, %zmm6{%k7}
	vaddpd        %zmm20, %zmm2, %zmm2
	vsubpd        %zmm22{cdab}, %zmm7, %zmm7{%k7}
	vaddpd        %zmm23, %zmm5, %zmm5
	vsubpd        %zmm23{cdab}, %zmm8, %zmm8{%k7}
	vsubpd        %zmm18{cdab}, %zmm9, %zmm9{%k7}
	vaddpd        %zmm21{cdab}, %zmm6, %zmm6{%k6}
	vsubpd        %zmm19{cdab}, %zmm10, %zmm10{%k7}
	vmovapd        576(%r10), %zmm21
	vaddpd        %zmm22{cdab}, %zmm7, %zmm7{%k6}
	vsubpd        %zmm20{cdab}, %zmm11, %zmm11{%k7}
	vmovapd        640(%r10), %zmm22
	vaddpd        %zmm23{cdab}, %zmm8, %zmm8{%k6}
	vaddpd        %zmm18{cdab}, %zmm9, %zmm9{%k6}
	vmovapd        704(%r10), %zmm23
	vpandd         272(%rsp){1to16} , %zmm21, %zmm21
	vmovapd        384(%r10), %zmm18
	vaddpd        %zmm19{cdab}, %zmm10, %zmm10{%k6}
	vpandd         272(%rsp){1to16} , %zmm22, %zmm22
	vmovapd        448(%r10), %zmm19
	vaddpd        %zmm20{cdab}, %zmm11, %zmm11{%k6}
	vpandd         272(%rsp){1to16} , %zmm18, %zmm18
	vmovapd        512(%r10), %zmm20
	vpandd         272(%rsp){1to16} , %zmm19, %zmm19
	movq          %r12, %r10
	vpandd         272(%rsp){1to16} , %zmm23, %zmm23
	movq           160(%rsp), %r12
	vpandd         272(%rsp){1to16} , %zmm20, %zmm20
	addq          %r15, %r12
	vsubpd        %zmm21, %zmm12, %zmm12
	movq           232(%rsp), %r15
	vsubpd        %zmm22, %zmm13, %zmm13
	vprefetch0       0(%r12)
	vsubpd        %zmm23, %zmm14, %zmm14
	vprefetch1     768(%r12)
	vaddpd        %zmm18, %zmm15, %zmm15
	vprefetch0      64(%r12)
	vaddpd        %zmm19, %zmm16, %zmm16
	vprefetch1     832(%r12)
	vaddpd        %zmm20, %zmm17, %zmm17
	vprefetch0     128(%r12)
	vmulpd        %zmm12, %zmm24, %zmm18
	vprefetch1     896(%r12)
	vmulpd        %zmm12{cdab}, %zmm24, %zmm28
	vprefetch0     192(%r12)
	vmulpd        %zmm15, %zmm24, %zmm21
	vprefetch1     960(%r12)
	vmulpd        %zmm15{cdab}, %zmm24, %zmm29
	vprefetch0     256(%r12)
	vfmadd231pd   %zmm13, %zmm25, %zmm18
	vmovapd        192(%r11), %zmm24
	vfmadd231pd        %zmm13{cdab}, %zmm25, %zmm28
	vprefetch1    1024(%r12)
	vfmadd231pd   %zmm16, %zmm25, %zmm21
	vprefetch0     320(%r12)
	vfmadd231pd        %zmm16{cdab}, %zmm25, %zmm29
	vprefetch1    1088(%r12)
	vfmadd231pd   %zmm14, %zmm26, %zmm18
	vmovapd        256(%r11), %zmm25
	vfmadd231pd        %zmm14{cdab}, %zmm26, %zmm28
	vprefetch0     384(%r12)
	vfmadd231pd   %zmm17, %zmm26, %zmm21
	vprefetch1    1152(%r12)
	vfmadd231pd        %zmm17{cdab}, %zmm26, %zmm29
	vprefetch0     448(%r12)
	vsubpd        %zmm18{cdab}, %zmm18, %zmm18{%k7}
	vmovapd        320(%r11), %zmm26
	vmulpd        %zmm12, %zmm24, %zmm19
	vprefetch1    1216(%r12)
	vsubpd        %zmm21{cdab}, %zmm21, %zmm21{%k7}
	vprefetch0     512(%r12)
	vaddpd        %zmm28{cdab}, %zmm28, %zmm18{%k6}
	vprefetch1    1280(%r12)
	vmulpd        %zmm12{cdab}, %zmm24, %zmm28
	vprefetch0     576(%r12)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm21{%k6}
	vprefetch1    1344(%r12)
	vmulpd        %zmm15, %zmm24, %zmm22
	vprefetch0     640(%r12)
	vmulpd        %zmm15{cdab}, %zmm24, %zmm29
	vprefetch1    1408(%r12)
	vfmadd231pd   %zmm13, %zmm25, %zmm19
	vmovapd        384(%r11), %zmm24
	vfmadd231pd        %zmm13{cdab}, %zmm25, %zmm28
	vprefetch0     704(%r12)
	vfmadd231pd   %zmm16, %zmm25, %zmm22
	vprefetch1    1472(%r12)
	vfmadd231pd        %zmm16{cdab}, %zmm25, %zmm29
	orq           %r12, %r12
	vfmadd231pd   %zmm14, %zmm26, %zmm19
	vmovapd        448(%r11), %zmm25
	vfmadd231pd        %zmm14{cdab}, %zmm26, %zmm28
	imulq         %r14, %r15
	vfmadd231pd   %zmm17, %zmm26, %zmm22
	movq          %r9, %r14
	vfmadd231pd        %zmm17{cdab}, %zmm26, %zmm29
	subq          %rax, %r14
	vsubpd        %zmm19{cdab}, %zmm19, %zmm19{%k7}
	vmovapd        512(%r11), %zmm26
	vmulpd        %zmm12, %zmm24, %zmm20
	addq          $576, %r11
	vsubpd        %zmm22{cdab}, %zmm22, %zmm22{%k7}
	vaddpd        %zmm28{cdab}, %zmm28, %zmm19{%k6}
	vmulpd        %zmm12{cdab}, %zmm24, %zmm28
	vaddpd        %zmm29{cdab}, %zmm29, %zmm22{%k6}
	vmovapd          0(%r10), %zmm12
	vmulpd        %zmm15, %zmm24, %zmm23
	vmulpd        %zmm15{cdab}, %zmm24, %zmm29
	vfmadd231pd   %zmm13, %zmm25, %zmm20
	vmovapd        192(%r10), %zmm15
	vfmadd231pd        %zmm13{cdab}, %zmm25, %zmm28
	vmovapd          0(%r11), %zmm24
	vfmadd231pd   %zmm16, %zmm25, %zmm23
	vmovapd         64(%r10), %zmm13
	vfmadd231pd        %zmm16{cdab}, %zmm25, %zmm29
	vfmadd231pd   %zmm14, %zmm26, %zmm20
	vmovapd        256(%r10), %zmm16
	vfmadd231pd        %zmm14{cdab}, %zmm26, %zmm28
	vmovapd         64(%r11), %zmm25
	vfmadd231pd   %zmm17, %zmm26, %zmm23
	vmovapd        128(%r10), %zmm14
	vfmadd231pd        %zmm17{cdab}, %zmm26, %zmm29
	vsubpd        %zmm20{cdab}, %zmm20, %zmm20{%k7}
	vmovapd        320(%r10), %zmm17
	vaddpd        %zmm18, %zmm0, %zmm0
	vmovapd        128(%r11), %zmm26
	vsubpd        %zmm23{cdab}, %zmm23, %zmm23{%k7}
	vaddpd        %zmm28{cdab}, %zmm28, %zmm20{%k6}
	vaddpd        %zmm19, %zmm1, %zmm1
	vaddpd        %zmm29{cdab}, %zmm29, %zmm23{%k6}
	vaddpd        %zmm21, %zmm3, %zmm3
	vaddpd        %zmm22, %zmm4, %zmm4
	vaddpd        %zmm21, %zmm6, %zmm6
	vaddpd        %zmm20, %zmm2, %zmm2
	vmovapd        576(%r10), %zmm21
	vaddpd        %zmm22, %zmm7, %zmm7
	vaddpd        %zmm23, %zmm5, %zmm5
	vmovapd        640(%r10), %zmm22
	vaddpd        %zmm23, %zmm8, %zmm8
	vsubpd        %zmm18, %zmm9, %zmm9
	vmovapd        704(%r10), %zmm23
	vsubpd        %zmm19, %zmm10, %zmm10
	vmovapd        384(%r10), %zmm18
	vsubpd        %zmm20, %zmm11, %zmm11
	vmovapd        448(%r10), %zmm19
	vpandd         280(%rsp){1to16} , %zmm21, %zmm21
	vmovapd        512(%r10), %zmm20
	vpandd         280(%rsp){1to16} , %zmm18, %zmm18
	movq          %r12, %r10
	vpandd         280(%rsp){1to16} , %zmm19, %zmm19
	movq           168(%rsp), %r12
	vpandd         280(%rsp){1to16} , %zmm20, %zmm20
	addq          %r15, %r12
	vpandd         280(%rsp){1to16} , %zmm22, %zmm22
	movq           240(%rsp), %r15
	vpandd         280(%rsp){1to16} , %zmm23, %zmm23
	vprefetch0       0(%r12)
	vaddpd        %zmm21, %zmm12, %zmm12
	vprefetch1     768(%r12)
	vsubpd        %zmm18, %zmm15, %zmm15
	vprefetch0      64(%r12)
	vaddpd        %zmm22, %zmm13, %zmm13
	vprefetch1     832(%r12)
	vaddpd        %zmm23, %zmm14, %zmm14
	vprefetch0     128(%r12)
	vsubpd        %zmm19, %zmm16, %zmm16
	vprefetch1     896(%r12)
	vsubpd        %zmm20, %zmm17, %zmm17
	vprefetch0     192(%r12)
	vmulpd        %zmm12, %zmm24, %zmm18
	vprefetch1     960(%r12)
	vmulpd        %zmm12{cdab}, %zmm24, %zmm28
	vprefetch0     256(%r12)
	vmulpd        %zmm15, %zmm24, %zmm21
	vprefetch1    1024(%r12)
	vmulpd        %zmm15{cdab}, %zmm24, %zmm29
	vprefetch0     320(%r12)
	vfmadd231pd   %zmm13, %zmm25, %zmm18
	vmovapd        192(%r11), %zmm24
	vfmadd231pd        %zmm13{cdab}, %zmm25, %zmm28
	vprefetch1    1088(%r12)
	vfmadd231pd   %zmm16, %zmm25, %zmm21
	vprefetch0     384(%r12)
	vfmadd231pd        %zmm16{cdab}, %zmm25, %zmm29
	vprefetch1    1152(%r12)
	vfmadd231pd   %zmm14, %zmm26, %zmm18
	vmovapd        256(%r11), %zmm25
	vfmadd231pd        %zmm14{cdab}, %zmm26, %zmm28
	vprefetch0     448(%r12)
	vfmadd231pd   %zmm17, %zmm26, %zmm21
	vprefetch1    1216(%r12)
	vfmadd231pd        %zmm17{cdab}, %zmm26, %zmm29
	vprefetch0     512(%r12)
	vsubpd        %zmm18{cdab}, %zmm18, %zmm18{%k7}
	vmovapd        320(%r11), %zmm26
	vmulpd        %zmm12, %zmm24, %zmm19
	vprefetch1    1280(%r12)
	vsubpd        %zmm21{cdab}, %zmm21, %zmm21{%k7}
	vprefetch0     576(%r12)
	vaddpd        %zmm28{cdab}, %zmm28, %zmm18{%k6}
	vprefetch1    1344(%r12)
	vmulpd        %zmm12{cdab}, %zmm24, %zmm28
	vprefetch0     640(%r12)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm21{%k6}
	vprefetch1    1408(%r12)
	vmulpd        %zmm15, %zmm24, %zmm22
	vprefetch0     704(%r12)
	vmulpd        %zmm15{cdab}, %zmm24, %zmm29
	vprefetch1    1472(%r12)
	vfmadd231pd   %zmm13, %zmm25, %zmm19
	vmovapd        384(%r11), %zmm24
	vfmadd231pd        %zmm13{cdab}, %zmm25, %zmm28
	orq           %r12, %r12
	vfmadd231pd   %zmm16, %zmm25, %zmm22
	imulq         %r14, %r15
	vfmadd231pd        %zmm16{cdab}, %zmm25, %zmm29
	movq          %r9, %r14
	vfmadd231pd   %zmm14, %zmm26, %zmm19
	vmovapd        448(%r11), %zmm25
	vfmadd231pd        %zmm14{cdab}, %zmm26, %zmm28
	subq          %rax, %r14
	vfmadd231pd   %zmm17, %zmm26, %zmm22
	vfmadd231pd        %zmm17{cdab}, %zmm26, %zmm29
	vsubpd        %zmm19{cdab}, %zmm19, %zmm19{%k7}
	vmovapd        512(%r11), %zmm26
	vmulpd        %zmm12, %zmm24, %zmm20
	addq          $576, %r11
	vsubpd        %zmm22{cdab}, %zmm22, %zmm22{%k7}
	vaddpd        %zmm28{cdab}, %zmm28, %zmm19{%k6}
	vmulpd        %zmm12{cdab}, %zmm24, %zmm28
	vaddpd        %zmm29{cdab}, %zmm29, %zmm22{%k6}
	vmovapd          0(%r10), %zmm12
	vmulpd        %zmm15, %zmm24, %zmm23
	vmulpd        %zmm15{cdab}, %zmm24, %zmm29
	vfmadd231pd   %zmm13, %zmm25, %zmm20
	vmovapd        192(%r10), %zmm15
	vfmadd231pd        %zmm13{cdab}, %zmm25, %zmm28
	vmovapd          0(%r11), %zmm24
	vfmadd231pd   %zmm16, %zmm25, %zmm23
	vmovapd         64(%r10), %zmm13
	vfmadd231pd        %zmm16{cdab}, %zmm25, %zmm29
	vfmadd231pd   %zmm14, %zmm26, %zmm20
	vmovapd        256(%r10), %zmm16
	vfmadd231pd        %zmm14{cdab}, %zmm26, %zmm28
	vmovapd         64(%r11), %zmm25
	vfmadd231pd   %zmm17, %zmm26, %zmm23
	vmovapd        128(%r10), %zmm14
	vfmadd231pd        %zmm17{cdab}, %zmm26, %zmm29
	vsubpd        %zmm20{cdab}, %zmm20, %zmm20{%k7}
	vmovapd        320(%r10), %zmm17
	vaddpd        %zmm18, %zmm0, %zmm0
	vmovapd        128(%r11), %zmm26
	vsubpd        %zmm23{cdab}, %zmm23, %zmm23{%k7}
	vaddpd        %zmm28{cdab}, %zmm28, %zmm20{%k6}
	vaddpd        %zmm19, %zmm1, %zmm1
	vaddpd        %zmm29{cdab}, %zmm29, %zmm23{%k6}
	vaddpd        %zmm21, %zmm3, %zmm3
	vaddpd        %zmm22, %zmm4, %zmm4
	vsubpd        %zmm21, %zmm6, %zmm6
	vaddpd        %zmm20, %zmm2, %zmm2
	vmovapd        576(%r10), %zmm21
	vsubpd        %zmm22, %zmm7, %zmm7
	vaddpd        %zmm23, %zmm5, %zmm5
	vmovapd        640(%r10), %zmm22
	vsubpd        %zmm23, %zmm8, %zmm8
	vaddpd        %zmm18, %zmm9, %zmm9
	vmovapd        704(%r10), %zmm23
	vaddpd        %zmm19, %zmm10, %zmm10
	vmovapd        384(%r10), %zmm18
	vaddpd        %zmm20, %zmm11, %zmm11
	vmovapd        448(%r10), %zmm19
	vpandd         288(%rsp){1to16} , %zmm21, %zmm21
	vmovapd        512(%r10), %zmm20
	vpandd         288(%rsp){1to16} , %zmm18, %zmm18
	movq          %r12, %r10
	vpandd         288(%rsp){1to16} , %zmm19, %zmm19
	movq           176(%rsp), %r12
	vpandd         288(%rsp){1to16} , %zmm20, %zmm20
	addq          %r15, %r12
	vpandd         288(%rsp){1to16} , %zmm22, %zmm22
	movq           248(%rsp), %r15
	vpandd         288(%rsp){1to16} , %zmm23, %zmm23
	vprefetch0       0(%r12)
	vsubpd        %zmm18{cdab}, %zmm12, %zmm12{%k7}
	vprefetch1     768(%r12)
	vsubpd        %zmm19{cdab}, %zmm13, %zmm13{%k7}
	vprefetch0      64(%r12)
	vsubpd        %zmm20{cdab}, %zmm14, %zmm14{%k7}
	vprefetch1     832(%r12)
	vaddpd        %zmm21{cdab}, %zmm15, %zmm15{%k7}
	vprefetch0     128(%r12)
	vaddpd        %zmm22{cdab}, %zmm16, %zmm16{%k7}
	vprefetch1     896(%r12)
	vaddpd        %zmm23{cdab}, %zmm17, %zmm17{%k7}
	vprefetch0     192(%r12)
	vaddpd        %zmm18{cdab}, %zmm12, %zmm12{%k6}
	vprefetch1     960(%r12)
	vaddpd        %zmm19{cdab}, %zmm13, %zmm13{%k6}
	vprefetch0     256(%r12)
	vaddpd        %zmm20{cdab}, %zmm14, %zmm14{%k6}
	vprefetch1    1024(%r12)
	vsubpd        %zmm21{cdab}, %zmm15, %zmm15{%k6}
	vprefetch0     320(%r12)
	vsubpd        %zmm22{cdab}, %zmm16, %zmm16{%k6}
	vprefetch1    1088(%r12)
	vsubpd        %zmm23{cdab}, %zmm17, %zmm17{%k6}
	vprefetch0     384(%r12)
	vmulpd        %zmm12, %zmm24, %zmm18
	vprefetch1    1152(%r12)
	vmulpd        %zmm12{cdab}, %zmm24, %zmm28
	vprefetch0     448(%r12)
	vfmadd231pd   %zmm13, %zmm25, %zmm18
	vprefetch1    1216(%r12)
	vmulpd        %zmm15, %zmm24, %zmm21
	vprefetch0     512(%r12)
	vmulpd        %zmm15{cdab}, %zmm24, %zmm29
	vprefetch1    1280(%r12)
	vfmadd231pd        %zmm13{cdab}, %zmm25, %zmm28
	vmovapd        192(%r11), %zmm24
	vfmadd231pd   %zmm16, %zmm25, %zmm21
	vprefetch0     576(%r12)
	vfmadd231pd        %zmm16{cdab}, %zmm25, %zmm29
	vprefetch1    1344(%r12)
	vfmadd231pd   %zmm14, %zmm26, %zmm18
	vmovapd        256(%r11), %zmm25
	vfmadd231pd        %zmm14{cdab}, %zmm26, %zmm28
	vprefetch0     640(%r12)
	vfmadd231pd   %zmm17, %zmm26, %zmm21
	vprefetch1    1408(%r12)
	vfmadd231pd        %zmm17{cdab}, %zmm26, %zmm29
	vprefetch0     704(%r12)
	vsubpd        %zmm18{cdab}, %zmm18, %zmm18{%k7}
	vmovapd        320(%r11), %zmm26
	vmulpd        %zmm12, %zmm24, %zmm19
	vprefetch1    1472(%r12)
	vsubpd        %zmm21{cdab}, %zmm21, %zmm21{%k7}
	orq           %r12, %r12
	vaddpd        %zmm28{cdab}, %zmm28, %zmm18{%k6}
	imulq         %r14, %r15
	vmulpd        %zmm12{cdab}, %zmm24, %zmm28
	movq          %r9, %r14
	vaddpd        %zmm29{cdab}, %zmm29, %zmm21{%k6}
	subq          %rax, %r14
	vmulpd        %zmm15, %zmm24, %zmm22
	addq          $1, %r14
	vmulpd        %zmm15{cdab}, %zmm24, %zmm29
	addq          $-1, %rax
	vfmadd231pd   %zmm13, %zmm25, %zmm19
	vmovapd        384(%r11), %zmm24
	vfmadd231pd        %zmm13{cdab}, %zmm25, %zmm28
	vfmadd231pd   %zmm16, %zmm25, %zmm22
	vfmadd231pd        %zmm16{cdab}, %zmm25, %zmm29
	vfmadd231pd   %zmm14, %zmm26, %zmm19
	vmovapd        448(%r11), %zmm25
	vfmadd231pd        %zmm14{cdab}, %zmm26, %zmm28
	vfmadd231pd   %zmm17, %zmm26, %zmm22
	vfmadd231pd        %zmm17{cdab}, %zmm26, %zmm29
	vsubpd        %zmm19{cdab}, %zmm19, %zmm19{%k7}
	vmovapd        512(%r11), %zmm26
	vmulpd        %zmm12, %zmm24, %zmm20
	addq          $576, %r11
	vsubpd        %zmm22{cdab}, %zmm22, %zmm22{%k7}
	vaddpd        %zmm28{cdab}, %zmm28, %zmm19{%k6}
	vmulpd        %zmm12{cdab}, %zmm24, %zmm28
	vaddpd        %zmm29{cdab}, %zmm29, %zmm22{%k6}
	vmovapd          0(%r10), %zmm12
	vmulpd        %zmm15, %zmm24, %zmm23
	vmulpd        %zmm15{cdab}, %zmm24, %zmm29
	vfmadd231pd   %zmm13, %zmm25, %zmm20
	vmovapd        192(%r10), %zmm15
	vfmadd231pd        %zmm13{cdab}, %zmm25, %zmm28
	vmovapd          0(%r11), %zmm24
	vfmadd231pd   %zmm16, %zmm25, %zmm23
	vmovapd         64(%r10), %zmm13
	vfmadd231pd        %zmm16{cdab}, %zmm25, %zmm29
	vfmadd231pd   %zmm14, %zmm26, %zmm20
	vmovapd        256(%r10), %zmm16
	vfmadd231pd        %zmm14{cdab}, %zmm26, %zmm28
	vmovapd         64(%r11), %zmm25
	vfmadd231pd   %zmm17, %zmm26, %zmm23
	vmovapd        128(%r10), %zmm14
	vfmadd231pd        %zmm17{cdab}, %zmm26, %zmm29
	vsubpd        %zmm20{cdab}, %zmm20, %zmm20{%k7}
	vmovapd        320(%r10), %zmm17
	vaddpd        %zmm18, %zmm0, %zmm0
	vmovapd        128(%r11), %zmm26
	vsubpd        %zmm23{cdab}, %zmm23, %zmm23{%k7}
	vaddpd        %zmm28{cdab}, %zmm28, %zmm20{%k6}
	vaddpd        %zmm19, %zmm1, %zmm1
	vaddpd        %zmm29{cdab}, %zmm29, %zmm23{%k6}
	vaddpd        %zmm21, %zmm3, %zmm3
	vaddpd        %zmm22, %zmm4, %zmm4
	vaddpd        %zmm18{cdab}, %zmm6, %zmm6{%k7}
	vaddpd        %zmm20, %zmm2, %zmm2
	vaddpd        %zmm19{cdab}, %zmm7, %zmm7{%k7}
	vaddpd        %zmm23, %zmm5, %zmm5
	vaddpd        %zmm20{cdab}, %zmm8, %zmm8{%k7}
	vsubpd        %zmm21{cdab}, %zmm9, %zmm9{%k7}
	vsubpd        %zmm18{cdab}, %zmm6, %zmm6{%k6}
	vsubpd        %zmm22{cdab}, %zmm10, %zmm10{%k7}
	vmovapd        384(%r10), %zmm18
	vsubpd        %zmm19{cdab}, %zmm7, %zmm7{%k6}
	vsubpd        %zmm23{cdab}, %zmm11, %zmm11{%k7}
	vmovapd        448(%r10), %zmm19
	vsubpd        %zmm20{cdab}, %zmm8, %zmm8{%k6}
	vaddpd        %zmm21{cdab}, %zmm9, %zmm9{%k6}
	vmovapd        512(%r10), %zmm20
	vpandd         296(%rsp){1to16} , %zmm18, %zmm18
	vmovapd        576(%r10), %zmm21
	vaddpd        %zmm22{cdab}, %zmm10, %zmm10{%k6}
	vpandd         296(%rsp){1to16} , %zmm19, %zmm19
	vmovapd        640(%r10), %zmm22
	vaddpd        %zmm23{cdab}, %zmm11, %zmm11{%k6}
	vpandd         296(%rsp){1to16} , %zmm20, %zmm20
	vmovapd        704(%r10), %zmm23
	vpandd         296(%rsp){1to16} , %zmm21, %zmm21
	movq          %r12, %r10
	vpandd         296(%rsp){1to16} , %zmm22, %zmm22
	movq           184(%rsp), %r12
	vpandd         296(%rsp){1to16} , %zmm23, %zmm23
	addq          %r15, %r12
	vaddpd        %zmm18{cdab}, %zmm12, %zmm12{%k7}
	movq           192(%rsp), %r15
	vaddpd        %zmm19{cdab}, %zmm13, %zmm13{%k7}
	vprefetch0       0(%r12)
	vaddpd        %zmm20{cdab}, %zmm14, %zmm14{%k7}
	vprefetch1     768(%r12)
	vsubpd        %zmm21{cdab}, %zmm15, %zmm15{%k7}
	vprefetch0      64(%r12)
	vsubpd        %zmm22{cdab}, %zmm16, %zmm16{%k7}
	vprefetch1     832(%r12)
	vsubpd        %zmm23{cdab}, %zmm17, %zmm17{%k7}
	vprefetch0     128(%r12)
	vsubpd        %zmm18{cdab}, %zmm12, %zmm12{%k6}
	vprefetch1     896(%r12)
	vsubpd        %zmm19{cdab}, %zmm13, %zmm13{%k6}
	vprefetch0     192(%r12)
	vsubpd        %zmm20{cdab}, %zmm14, %zmm14{%k6}
	vprefetch1     960(%r12)
	vaddpd        %zmm21{cdab}, %zmm15, %zmm15{%k6}
	vprefetch0     256(%r12)
	vaddpd        %zmm22{cdab}, %zmm16, %zmm16{%k6}
	vprefetch1    1024(%r12)
	vaddpd        %zmm23{cdab}, %zmm17, %zmm17{%k6}
	vprefetch0     320(%r12)
	vmulpd        %zmm12, %zmm24, %zmm18
	vprefetch1    1088(%r12)
	vmulpd        %zmm12{cdab}, %zmm24, %zmm28
	vprefetch0     384(%r12)
	vfmadd231pd   %zmm13, %zmm25, %zmm18
	vprefetch1    1152(%r12)
	vmulpd        %zmm15, %zmm24, %zmm21
	vprefetch0     448(%r12)
	vmulpd        %zmm15{cdab}, %zmm24, %zmm29
	vprefetch1    1216(%r12)
	vfmadd231pd        %zmm13{cdab}, %zmm25, %zmm28
	vmovapd        192(%r11), %zmm24
	vfmadd231pd   %zmm16, %zmm25, %zmm21
	vprefetch0     512(%r12)
	vfmadd231pd        %zmm16{cdab}, %zmm25, %zmm29
	vprefetch1    1280(%r12)
	vfmadd231pd   %zmm14, %zmm26, %zmm18
	vmovapd        256(%r11), %zmm25
	vfmadd231pd        %zmm14{cdab}, %zmm26, %zmm28
	vprefetch0     576(%r12)
	vfmadd231pd   %zmm17, %zmm26, %zmm21
	vprefetch1    1344(%r12)
	vfmadd231pd        %zmm17{cdab}, %zmm26, %zmm29
	vprefetch0     640(%r12)
	vsubpd        %zmm18{cdab}, %zmm18, %zmm18{%k7}
	vmovapd        320(%r11), %zmm26
	vmulpd        %zmm12, %zmm24, %zmm19
	vprefetch1    1408(%r12)
	vsubpd        %zmm21{cdab}, %zmm21, %zmm21{%k7}
	vprefetch0     704(%r12)
	vaddpd        %zmm28{cdab}, %zmm28, %zmm18{%k6}
	vprefetch1    1472(%r12)
	vmulpd        %zmm12{cdab}, %zmm24, %zmm28
	orq           %r12, %r12
	vaddpd        %zmm29{cdab}, %zmm29, %zmm21{%k6}
	imulq         %r14, %r15
	vmulpd        %zmm15, %zmm24, %zmm22
	vmulpd        %zmm15{cdab}, %zmm24, %zmm29
	vfmadd231pd   %zmm13, %zmm25, %zmm19
	vmovapd        384(%r11), %zmm24
	vfmadd231pd        %zmm13{cdab}, %zmm25, %zmm28
	vfmadd231pd   %zmm16, %zmm25, %zmm22
	vfmadd231pd        %zmm16{cdab}, %zmm25, %zmm29
	vfmadd231pd   %zmm14, %zmm26, %zmm19
	vmovapd        448(%r11), %zmm25
	vfmadd231pd        %zmm14{cdab}, %zmm26, %zmm28
	vfmadd231pd   %zmm17, %zmm26, %zmm22
	vfmadd231pd        %zmm17{cdab}, %zmm26, %zmm29
	vsubpd        %zmm19{cdab}, %zmm19, %zmm19{%k7}
	vmovapd        512(%r11), %zmm26
	vmulpd        %zmm12, %zmm24, %zmm20
	addq          $576, %r11
	vsubpd        %zmm22{cdab}, %zmm22, %zmm22{%k7}
	vaddpd        %zmm28{cdab}, %zmm28, %zmm19{%k6}
	vmulpd        %zmm12{cdab}, %zmm24, %zmm28
	vaddpd        %zmm29{cdab}, %zmm29, %zmm22{%k6}
	vmovapd          0(%r10), %zmm12
	vmulpd        %zmm15, %zmm24, %zmm23
	vmulpd        %zmm15{cdab}, %zmm24, %zmm29
	vfmadd231pd   %zmm13, %zmm25, %zmm20
	vmovapd        192(%r10), %zmm15
	vfmadd231pd        %zmm13{cdab}, %zmm25, %zmm28
	vmovapd          0(%r11), %zmm24
	vfmadd231pd   %zmm16, %zmm25, %zmm23
	vmovapd         64(%r10), %zmm13
	vfmadd231pd        %zmm16{cdab}, %zmm25, %zmm29
	vfmadd231pd   %zmm14, %zmm26, %zmm20
	vmovapd        256(%r10), %zmm16
	vfmadd231pd        %zmm14{cdab}, %zmm26, %zmm28
	vmovapd         64(%r11), %zmm25
	vfmadd231pd   %zmm17, %zmm26, %zmm23
	vmovapd        128(%r10), %zmm14
	vfmadd231pd        %zmm17{cdab}, %zmm26, %zmm29
	vsubpd        %zmm20{cdab}, %zmm20, %zmm20{%k7}
	vmovapd        320(%r10), %zmm17
	vaddpd        %zmm18, %zmm0, %zmm0
	vmovapd        128(%r11), %zmm26
	vsubpd        %zmm23{cdab}, %zmm23, %zmm23{%k7}
	vaddpd        %zmm28{cdab}, %zmm28, %zmm20{%k6}
	vaddpd        %zmm19, %zmm1, %zmm1
	vaddpd        %zmm29{cdab}, %zmm29, %zmm23{%k6}
	vaddpd        %zmm21, %zmm3, %zmm3
	vaddpd        %zmm22, %zmm4, %zmm4
	vsubpd        %zmm18{cdab}, %zmm6, %zmm6{%k7}
	vaddpd        %zmm20, %zmm2, %zmm2
	vsubpd        %zmm19{cdab}, %zmm7, %zmm7{%k7}
	vaddpd        %zmm23, %zmm5, %zmm5
	vsubpd        %zmm20{cdab}, %zmm8, %zmm8{%k7}
	vaddpd        %zmm21{cdab}, %zmm9, %zmm9{%k7}
	vaddpd        %zmm18{cdab}, %zmm6, %zmm6{%k6}
	vaddpd        %zmm22{cdab}, %zmm10, %zmm10{%k7}
	vmovapd        384(%r10), %zmm18
	vaddpd        %zmm19{cdab}, %zmm7, %zmm7{%k6}
	vaddpd        %zmm23{cdab}, %zmm11, %zmm11{%k7}
	vmovapd        448(%r10), %zmm19
	vaddpd        %zmm20{cdab}, %zmm8, %zmm8{%k6}
	vsubpd        %zmm21{cdab}, %zmm9, %zmm9{%k6}
	vmovapd        512(%r10), %zmm20
	vpandd         304(%rsp){1to16} , %zmm18, %zmm18
	vmovapd        576(%r10), %zmm21
	vsubpd        %zmm22{cdab}, %zmm10, %zmm10{%k6}
	vpandd         304(%rsp){1to16} , %zmm19, %zmm19
	vmovapd        640(%r10), %zmm22
	vsubpd        %zmm23{cdab}, %zmm11, %zmm11{%k6}
	vpandd         304(%rsp){1to16} , %zmm20, %zmm20
	vmovapd        704(%r10), %zmm23
	vpandd         304(%rsp){1to16} , %zmm21, %zmm21
	movq          %r12, %r10
	vpandd         304(%rsp){1to16} , %zmm22, %zmm22
	movq           128(%rsp), %r12
	vpandd         304(%rsp){1to16} , %zmm23, %zmm23
	addq          %r15, %r12
	vaddpd        %zmm18, %zmm12, %zmm12
	movq            80(%rsp), %r15
	vaddpd        %zmm19, %zmm13, %zmm13
	vprefetch0       0(%r12)
	vaddpd        %zmm20, %zmm14, %zmm14
	vprefetch1     768(%r12)
	vaddpd        %zmm21, %zmm15, %zmm15
	vprefetch0      64(%r12)
	vaddpd        %zmm22, %zmm16, %zmm16
	vprefetch1     832(%r12)
	vaddpd        %zmm23, %zmm17, %zmm17
	vprefetch0     128(%r12)
	vmulpd        %zmm12, %zmm24, %zmm18
	vprefetch1     896(%r12)
	vmulpd        %zmm12{cdab}, %zmm24, %zmm28
	vprefetch0     192(%r12)
	vmulpd        %zmm15, %zmm24, %zmm21
	vprefetch1     960(%r12)
	vmulpd        %zmm15{cdab}, %zmm24, %zmm29
	vprefetch0     256(%r12)
	vfmadd231pd   %zmm13, %zmm25, %zmm18
	vmovapd        192(%r11), %zmm24
	vfmadd231pd        %zmm13{cdab}, %zmm25, %zmm28
	vprefetch1    1024(%r12)
	vfmadd231pd   %zmm16, %zmm25, %zmm21
	vprefetch0     320(%r12)
	vfmadd231pd        %zmm16{cdab}, %zmm25, %zmm29
	vprefetch1    1088(%r12)
	vfmadd231pd   %zmm14, %zmm26, %zmm18
	vmovapd        256(%r11), %zmm25
	vfmadd231pd        %zmm14{cdab}, %zmm26, %zmm28
	vprefetch0     384(%r12)
	vfmadd231pd   %zmm17, %zmm26, %zmm21
	vprefetch1    1152(%r12)
	vfmadd231pd        %zmm17{cdab}, %zmm26, %zmm29
	vprefetch0     448(%r12)
	vsubpd        %zmm18{cdab}, %zmm18, %zmm18{%k7}
	vmovapd        320(%r11), %zmm26
	vmulpd        %zmm12, %zmm24, %zmm19
	vprefetch1    1216(%r12)
	vsubpd        %zmm21{cdab}, %zmm21, %zmm21{%k7}
	vprefetch0     512(%r12)
	vaddpd        %zmm28{cdab}, %zmm28, %zmm18{%k6}
	vprefetch1    1280(%r12)
	vmulpd        %zmm12{cdab}, %zmm24, %zmm28
	vprefetch0     576(%r12)
	vaddpd        %zmm29{cdab}, %zmm29, %zmm21{%k6}
	vprefetch1    1344(%r12)
	vmulpd        %zmm15, %zmm24, %zmm22
	vprefetch0     640(%r12)
	vmulpd        %zmm15{cdab}, %zmm24, %zmm29
	vprefetch1    1408(%r12)
	vfmadd231pd   %zmm13, %zmm25, %zmm19
	vmovapd        384(%r11), %zmm24
	vfmadd231pd        %zmm13{cdab}, %zmm25, %zmm28
	vprefetch0     704(%r12)
	vfmadd231pd   %zmm16, %zmm25, %zmm22
	vprefetch1    1472(%r12)
	vfmadd231pd        %zmm16{cdab}, %zmm25, %zmm29
	vfmadd231pd   %zmm14, %zmm26, %zmm19
	vmovapd        448(%r11), %zmm25
	vfmadd231pd        %zmm14{cdab}, %zmm26, %zmm28
	vfmadd231pd   %zmm17, %zmm26, %zmm22
	vfmadd231pd        %zmm17{cdab}, %zmm26, %zmm29
	vsubpd        %zmm19{cdab}, %zmm19, %zmm19{%k7}
	vmovapd        512(%r11), %zmm26
	vmulpd        %zmm12, %zmm24, %zmm20
	addq          $576, %r11
	vsubpd        %zmm22{cdab}, %zmm22, %zmm22{%k7}
	vaddpd        %zmm28{cdab}, %zmm28, %zmm19{%k6}
	vmulpd        %zmm12{cdab}, %zmm24, %zmm28
	vaddpd        %zmm29{cdab}, %zmm29, %zmm22{%k6}
	vmovapd          0(%r10), %zmm12
	vmulpd        %zmm15, %zmm24, %zmm23
	vmulpd        %zmm15{cdab}, %zmm24, %zmm29
	vfmadd231pd   %zmm13, %zmm25, %zmm20
	vmovapd        192(%r10), %zmm15
	vfmadd231pd        %zmm13{cdab}, %zmm25, %zmm28
	vmovapd          0(%r11), %zmm24
	vfmadd231pd   %zmm16, %zmm25, %zmm23
	vmovapd         64(%r10), %zmm13
	vfmadd231pd        %zmm16{cdab}, %zmm25, %zmm29
	vfmadd231pd   %zmm14, %zmm26, %zmm20
	vmovapd        256(%r10), %zmm16
	vfmadd231pd        %zmm14{cdab}, %zmm26, %zmm28
	vmovapd         64(%r11), %zmm25
	vfmadd231pd   %zmm17, %zmm26, %zmm23
	vmovapd        128(%r10), %zmm14
	vfmadd231pd        %zmm17{cdab}, %zmm26, %zmm29
	vsubpd        %zmm20{cdab}, %zmm20, %zmm20{%k7}
	vmovapd        320(%r10), %zmm17
	vaddpd        %zmm18, %zmm0, %zmm0
	vmovapd        128(%r11), %zmm26
	vsubpd        %zmm23{cdab}, %zmm23, %zmm23{%k7}
	vaddpd        %zmm28{cdab}, %zmm28, %zmm20{%k6}
	vaddpd        %zmm19, %zmm1, %zmm1
	vaddpd        %zmm29{cdab}, %zmm29, %zmm23{%k6}
	vaddpd        %zmm21, %zmm3, %zmm3
	vaddpd        %zmm22, %zmm4, %zmm4
	vaddpd        %zmm18, %zmm6, %zmm6
	vaddpd        %zmm20, %zmm2, %zmm2
	vmovapd        384(%r10), %zmm18
	vaddpd        %zmm19, %zmm7, %zmm7
	vaddpd        %zmm23, %zmm5, %zmm5
	vmovapd        448(%r10), %zmm19
	vaddpd        %zmm20, %zmm8, %zmm8
	vaddpd        %zmm21, %zmm9, %zmm9
	vmovapd        512(%r10), %zmm20
	vaddpd        %zmm22, %zmm10, %zmm10
	vmovapd        576(%r10), %zmm21
	vaddpd        %zmm23, %zmm11, %zmm11
	vmovapd        640(%r10), %zmm22
	vpandd         312(%rsp){1to16} , %zmm18, %zmm18
	vmovapd        704(%r10), %zmm23
	vpandd         312(%rsp){1to16} , %zmm19, %zmm19
	vpandd         312(%rsp){1to16} , %zmm20, %zmm20
	vpandd         312(%rsp){1to16} , %zmm21, %zmm21
	vpandd         312(%rsp){1to16} , %zmm22, %zmm22
	vpandd         312(%rsp){1to16} , %zmm23, %zmm23
	vsubpd        %zmm18, %zmm12, %zmm12
	vsubpd        %zmm19, %zmm13, %zmm13
	vsubpd        %zmm20, %zmm14, %zmm14
	vsubpd        %zmm21, %zmm15, %zmm15
	vsubpd        %zmm22, %zmm16, %zmm16
	vsubpd        %zmm23, %zmm17, %zmm17
	vmulpd        %zmm12, %zmm24, %zmm18
	vmulpd        %zmm12{cdab}, %zmm24, %zmm28
	vmulpd        %zmm15, %zmm24, %zmm21
	vmulpd        %zmm15{cdab}, %zmm24, %zmm29
	vfmadd231pd   %zmm13, %zmm25, %zmm18
	vmovapd        192(%r11), %zmm24
	vfmadd231pd        %zmm13{cdab}, %zmm25, %zmm28
	vfmadd231pd   %zmm16, %zmm25, %zmm21
	vfmadd231pd        %zmm16{cdab}, %zmm25, %zmm29
	vfmadd231pd   %zmm14, %zmm26, %zmm18
	vmovapd        256(%r11), %zmm25
	vfmadd231pd        %zmm14{cdab}, %zmm26, %zmm28
	vfmadd231pd   %zmm17, %zmm26, %zmm21
	vfmadd231pd        %zmm17{cdab}, %zmm26, %zmm29
	vsubpd        %zmm18{cdab}, %zmm18, %zmm18{%k7}
	vmovapd        320(%r11), %zmm26
	vmulpd        %zmm12, %zmm24, %zmm19
	vsubpd        %zmm21{cdab}, %zmm21, %zmm21{%k7}
	vaddpd        %zmm28{cdab}, %zmm28, %zmm18{%k6}
	vmulpd        %zmm12{cdab}, %zmm24, %zmm28
	vaddpd        %zmm29{cdab}, %zmm29, %zmm21{%k6}
	vmulpd        %zmm15, %zmm24, %zmm22
	vmulpd        %zmm15{cdab}, %zmm24, %zmm29
	vfmadd231pd   %zmm13, %zmm25, %zmm19
	vmovapd        384(%r11), %zmm24
	vfmadd231pd        %zmm13{cdab}, %zmm25, %zmm28
	vfmadd231pd   %zmm16, %zmm25, %zmm22
	vfmadd231pd        %zmm16{cdab}, %zmm25, %zmm29
	vfmadd231pd   %zmm14, %zmm26, %zmm19
	vmovapd        448(%r11), %zmm25
	vfmadd231pd        %zmm14{cdab}, %zmm26, %zmm28
	vfmadd231pd   %zmm17, %zmm26, %zmm22
	vfmadd231pd        %zmm17{cdab}, %zmm26, %zmm29
	vsubpd        %zmm19{cdab}, %zmm19, %zmm19{%k7}
	vmovapd        512(%r11), %zmm26
	vmulpd        %zmm12, %zmm24, %zmm20
	addq          $576, %r11
	vsubpd        %zmm22{cdab}, %zmm22, %zmm22{%k7}
	vaddpd        %zmm28{cdab}, %zmm28, %zmm19{%k6}
	vmulpd        %zmm12{cdab}, %zmm24, %zmm28
	vaddpd        %zmm29{cdab}, %zmm29, %zmm22{%k6}
	vmulpd        %zmm15, %zmm24, %zmm23
	vmulpd        %zmm15{cdab}, %zmm24, %zmm29
	vfmadd231pd   %zmm13, %zmm25, %zmm20
	vfmadd231pd        %zmm13{cdab}, %zmm25, %zmm28
	vfmadd231pd   %zmm16, %zmm25, %zmm23
	vfmadd231pd        %zmm16{cdab}, %zmm25, %zmm29
	vfmadd231pd   %zmm14, %zmm26, %zmm20
	vfmadd231pd        %zmm14{cdab}, %zmm26, %zmm28
	vfmadd231pd   %zmm17, %zmm26, %zmm23
	vfmadd231pd        %zmm17{cdab}, %zmm26, %zmm29
	vsubpd        %zmm20{cdab}, %zmm20, %zmm20{%k7}
	vaddpd        %zmm18, %zmm0, %zmm0
	vsubpd        %zmm23{cdab}, %zmm23, %zmm23{%k7}
	vaddpd        %zmm28{cdab}, %zmm28, %zmm20{%k6}
	vaddpd        %zmm19, %zmm1, %zmm1
	vaddpd        %zmm29{cdab}, %zmm29, %zmm23{%k6}
	vmovnrngoapd  %zmm0,    0(%r15)
	vaddpd        %zmm21, %zmm3, %zmm3
	vaddpd        %zmm22, %zmm4, %zmm4
	vsubpd        %zmm18, %zmm6, %zmm6
	vmovnrngoapd  %zmm1,   64(%r15)
	vaddpd        %zmm20, %zmm2, %zmm2
	vsubpd        %zmm19, %zmm7, %zmm7
	vmovnrngoapd  %zmm3,  192(%r15)
	vaddpd        %zmm23, %zmm5, %zmm5
	vmovnrngoapd  %zmm4,  256(%r15)
	vsubpd        %zmm20, %zmm8, %zmm8
	vmovnrngoapd  %zmm6,  384(%r15)
	vsubpd        %zmm21, %zmm9, %zmm9
	vmovnrngoapd  %zmm2,  128(%r15)
	vsubpd        %zmm22, %zmm10, %zmm10
	vmovnrngoapd  %zmm7,  448(%r15)
	vsubpd        %zmm23, %zmm11, %zmm11
	vmovnrngoapd  %zmm5,  320(%r15)
	vmovnrngoapd  %zmm8,  512(%r15)
	vmovnrngoapd  %zmm9,  576(%r15)
	vmovnrngoapd  %zmm10,  640(%r15)
	vmovnrngoapd  %zmm11,  704(%r15)
	addq          $768, %r15
	orq           %rax, %rax
	movq          %r15,   80(%rsp)
	jg     _vmx_dwf_ext_lab26
	movq            64(%rsp), %r13
	addq          $-1, %r13
	orq           %r13, %r13
	jg     _vmx_dwf_ext_lab1
_vmx_dwf_ext_lab0:
	addq          $0, %rsp
	movq          2680(%rsp), %r15
	movq          2672(%rsp), %r14
	movq          2664(%rsp), %r13
	movq          2656(%rsp), %r12
	movq          2632(%rsp), %rbp
	movq          2624(%rsp), %rbx
	addq          $2688, %rsp

	# ----- Exit_Routine -----
	movq     (%rsp), %rsp
	popq     %r10
	ret
.align 16,0x90
.type vmx_dwf_ext,@function
.size vmx_dwf_ext,.-vmx_dwf_ext
.data
