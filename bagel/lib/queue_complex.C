/*
 *
 *  Copyright Peter Boyle and Glasgow University, 2000.
 *
 *  It is provided under the GNU pubic License V2
 *  It is provided as is and is not guaranteed fit for any purpose.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "processor.h"



/* {...COMPLEX queueing routines */

extern struct processor *PROC;
extern int pregs[];
extern int pimms[];

#define  MAX(a, b)      ((a) > (b) ? (a) : (b))
#define  MIN(a, b)      ((a) < (b) ? (a) : (b))

static int have_fmadd;
static int fp_latency;

struct rotating_reg *CMADregs = NULL;

int have_hummer(void)
{
	int have_hum = 0;
	int i;
	for (i =0;i < PROC->npipes; i++){
		if (PROC->pipes[i].issuegroups & SIMDHUMMERMASK) have_hum = 1;
	}
	return have_hum;
}

int have_permute(void)
{
	int i, have;

	have = 0;
	for (i = 0; i < PROC->npipes; i++)
	{
		if (PROC->pipes[i].issuegroups & SIMDPERMUTEMASK)
			have = 1;
	}
	return have;
}


/**
* \return Returns 1 if MIC is supported, returns 0 otherwise.
*/
int have_mic()
{
	int i, have;

	have = 0;
	for (i = 0; i < PROC->npipes; i++)
	{
		if ((PROC->pipes[i].issuegroups & SIMDMICMASK) || (PROC->pipes[i].issuegroups & SIMDMICSWIZMASK))
			have = 1;
	}
	return have;
}


int nsimd(void)
{
	return PROC->nsimd;
}



int reg_tmp = -1;
int i_addr = -1;
int which_is_loaded = -1;
int imm_0 = -1;
int imm_1 = -1;
int imm_2 = -1;


void complex_constants_prepare(int creg, int i_ea)  // XXX: MIC support?
{
	reg_tmp = creg;
	i_addr  = i_ea;

	if (imm_0 == -1) {
		imm_0 = def_offset(0,Byte,          "simd_zero_offset");
		imm_1 = def_offset(2*nsimd(),Double,"simd_one_offset");
		imm_2 = def_offset(4*nsimd(),Double,"simd_two_offset");
	}
	complex_load(reg_tmp,imm_0,i_addr,Double);
	which_is_loaded = cmplx_i;
}


void complex_constants_assert(int constant)
{
	if (which_is_loaded == constant)
		return;

	if (constant == cmplx_zero)
	{
		printf("complex_constants_assert: Not supported cmplx_zero");
	}
	fprintf(stderr, "complex_constants_assert : reg_tmp = %d\n", reg_tmp);
	exit(-1);
}


void setup_cmadds(void)
{
	int i;

	fp_latency = 0;
	have_fmadd = 0;
	test_proc();

	for (i = 0; i < PROC->npipes; i++)
	{
		if (PROC->pipes[i].issuegroups & FMACMASK)
		{
			have_fmadd = 1;
		}
		fp_latency = PROC->latency[FMULPIPE];
	}

	if (!have_fmadd)
	{
		CMADregs = create_rotating_reg(Fregs, MAX(1+fp_latency, 2), "CMADD");
	}
}


/**
 * Element-by-element SIMD MULTIPLY: e = a * c.
 * (Parameters are indexes in register files (see processor.C: make_inst()).)
 */
void simd_mul(int e, int a, int c)
{
	if (have_hummer())
	{
		make_inst(SIMDHUMMERPIPE, VMUL, e, a, c);
	}
	else if (have_mic())
	{
		make_inst(SIMDMICPIPE, MIC_MUL, e, a, c);
	}
	else
	{
		for (int v = 0; v < nsimd(); v++)
		{
			int f = e+1;  /*RE/IM innermost for hummer style SIMD*/
			int b = a+1;
			int d = c+1;
			queue_fmul(e, a, c);
			queue_fmul(f, b, d);
			e+=2; a+=2; c+=2;
		}
	}
}


/**
 * Element-by-element SIMD ADD: e = a + c
 */
void simd_add(int e, int a, int c)
{
	if (have_hummer())
	{
		make_inst(SIMDHUMMERPIPE, VADD, e, a, c);
	}
	else if (have_mic())
	{
		make_inst(SIMDMICPIPE, MIC_ADD, e, a, c);
	}
	else
	{
		for (int v = 0; v < nsimd(); v++)
		{
			int f = e+1;  /*RE/IM innermost for hummer style SIMD*/
			int b = a+1;
			int d = c+1;
			queue_fadd(e, a, c);
			queue_fadd(f, b, d);
			e+=2; a+=2; c+=2;
		}
	}
}


/**
 * Element-by-element SIMD multiply-add: e = (a*c) + g.
 */
void simd_madd(int e, int a, int c, int g)
{
	if (have_hummer())
	{
		make_inst(SIMDHUMMERPIPE, VMADD, e, a, c, g);
	}
	else if (have_mic())
	{
		// MIC's VFMADD has 3 operands only. Must reduce 4 => 3 ...
		if (e == a)      make_inst(SIMDMICPIPE, MIC_FMADD213, a, c, g);
		else if (e == c) make_inst(SIMDMICPIPE, MIC_FMADD213, c, a, g);
		else if (e == g) make_inst(SIMDMICPIPE, MIC_FMADD231, g, a, c);
		else
		{
			// NB: must copy g=>e to provide non-destructive 4-operand MADD
			make_inst(SIMDMICPIPE, MIC_MOV, e, g);  // copy: g => e
			make_inst(SIMDMICPIPE, MIC_FMADD231, e, a, c);
		}
	}
	else
	{
		for (int v = 0; v < nsimd(); v++)
		{
			int f = e+1;  /*RE/IM innermost for hummer style SIMD*/
			int b = a+1;
			int d = c+1;
			int h = g+1;
			queue_fmadd(e, a, c, g);
			queue_fmadd(f, b, d, h);
			e+=2; a+=2; c+=2; g+=2;
		}
	}
}


/**
 * e = a*c - g
 */
void simd_msub(int e, int a, int c, int g)
{
	if (have_mic())
	{
		// MIC's VFMSUB has 3 operands only. Must reduce 4 => 3 ...
		if (e == a)      make_inst(SIMDMICPIPE, MIC_FMSUB213, a, c, g);
		else if (e == c) make_inst(SIMDMICPIPE, MIC_FMSUB213, c, a, g);
		else if (e == g) make_inst(SIMDMICPIPE, MIC_FMSUB231, g, a, c);
		else
		{
			// NB: must copy g=>e to provide non-destructive 4-operand MADD
			make_inst(SIMDMICPIPE, MIC_MOV, e, g);  // copy: g => e
			make_inst(SIMDMICPIPE, MIC_FMSUB231, e, a, c);
		}
	}
	else {
		fprintf(stderr, "FIXME: %s not implemented for %s.\n", __FUNCTION__, procargs[PROC->id]);
		abort();
	}
}


/**
 * \brief Sets vector reg to zero.
 * \param[out] a Vector register set to zero.
 */
void simd_zero(int a)
{
	if (have_mic()) {
		make_inst(SIMDMICPIPE, MIC_XOR, a, a, a);
	}
	else {
		fprintf(stderr, "FIXME: %s not implemented for %s.\n", __FUNCTION__, procargs[PROC->id]);
		abort();
	}
}


/**
 * Horizontal ADD of \p v and store to offhandle(basereg) = v[0]+v[1]+...v[n-1]
 */
void simd_reduce_add(int basereg, int offhandle, int v, enum Datum type)
{
	if (have_mic()) {
		int t, s, pattern;

		t = allocate_tmp_reg(Cregs);
		s = allocate_tmp_reg(Cregs);
		pattern = def_offset(MIC_PERM_DCDC, Byte, "PERM_DCDC");
		make_inst(SIMDMICSWIZPIPE, MIC_ADD_SWIZ, t, v, v, MIC_SWIZ_CDAB);
		make_inst(SIMDMICSWIZPIPE, MIC_ADD_SWIZ, t, t, t, MIC_SWIZ_BADC);
		make_inst(SIMDPERMUTEPIPE, VPERMCOARSE, s, pattern, t);
		make_inst(SIMDMICPIPE, MIC_ADD, t, t, s);
		queue_fstore(t, offhandle, basereg, type);
		free_reg(Cregs, t);
		free_reg(Cregs, s);
	}
	else {
		fprintf(stderr, "FIXME: %s not implemented for %s.\n", __FUNCTION__, procargs[PROC->id]);
		abort();
	}
}


/**
 * Element-by-element NEGATE-MULTIPLY AND ADD: e = g - (a*c).
 */
void simd_nmadd(int e, int a, int c, int g)
{
	if (have_hummer())
	{
		make_inst(SIMDHUMMERPIPE, VNMSUB, e, a, c, g);
	}
	else if (have_mic())
	{
		if (e == a)      make_inst(SIMDMICPIPE, MIC_FNMADD213, a, c, g);
		else if (e == c) make_inst(SIMDMICPIPE, MIC_FNMADD213, c, a, g);
		else if (e == g) make_inst(SIMDMICPIPE, MIC_FNMADD231, g, a, c);
		else
		{
			make_inst(SIMDMICPIPE, MIC_MOV, e, g);
			make_inst(SIMDMICPIPE, MIC_FNMADD231, e, a, c);
		}
	}
	else
	{
		for (int v = 0; v < nsimd(); v++)
		{
			int f = e+1;  /*RE/IM innermost for hummer style SIMD*/
			int b = a+1;
			int d = c+1;
			int h = g+1;
			queue_fnmsub(e, a, c, g);
			queue_fnmsub(f, b, d, h);
			e+=2; a+=2; c+=2; g+=2;
		}
	}
}


/**
 * Queues vector complex MUL: (e + if) = (a+ib)*(c+id)
 */
void complex_mul(int e, int a, int c)
{
	if (have_hummer())
	{
		// complex_constants_assert(cmplx_zero);
		make_inst(SIMDHUMMERPIPE, VMUL_RR_RI, e, a, c, reg_tmp);
		make_inst(SIMDHUMMERPIPE, VMADD_MII_IR, e, a, c, e);
	}
	else if (have_mic())
	{
		int t;

		t = allocate_tmp_reg(Fregs);
		make_inst(SIMDMICPIPE, MIC_MUL, e, a, c);  // VMUL_RR_II
		make_inst(SIMDMICSWIZPIPE, MIC_MUL_SWIZ, t, a, c, MIC_SWIZ_CDAB);  // VMUL_RI_IR
		make_inst(SIMDMICSWIZPIPE, MIC_SUB_MASK_SWIZ, e, MIC_MASK_ODD,  e, e, MIC_SWIZ_CDAB);  // (RR-II,??)
		make_inst(SIMDMICSWIZPIPE, MIC_ADD_MASK_SWIZ, e, MIC_MASK_EVEN, t, t, MIC_SWIZ_CDAB);  // (RR-II,RI+IR)
		free_reg(Fregs, t);
	}
	else
	{
		for (int v = 0; v < nsimd(); v++)
		{
			int f = e+1;  /*RE/IM innermost for hummer style SIMD*/
			int b = a+1;
			int d = c+1;
			queue_cmul(e,f,a,b,c,d);
			e+=2; a+=2; c+=2;
		}
	}
}


// scalar complex mul: (e+if)=(a+ib)*(c+id)
void queue_cmul(int e, int f, int a, int b, int c, int d)
{
	if (have_fmadd)
	{
		if (have_mic()) {
			make_inst(FMULPIPE, FMUL, e, a, c);
			make_inst(FMULPIPE, FMUL, f, a, d);
			make_inst(FMACPIPE, FNMADD231, e, b, d);  // e -= b*d
			make_inst(FMACPIPE, FMADD231,  f, b, c);  // f += b*c
		}
		else {
			make_inst(FMULPIPE,FMUL,e,a,c);
			make_inst(FMULPIPE,FMUL,f,a,d);
			make_inst(FMACPIPE,FNMSUB,e,b,d,e);
			make_inst(FMACPIPE,FMADD,f,b,c,f);
		}
	}
	else
	{
		int t;

		make_inst(FMULPIPE,FMUL,e,a,c);
		make_inst(FMULPIPE,FMUL,f,a,d);

		t = get_rotating_register(CMADregs);
		make_inst(FMULPIPE,FMUL,t,b,d);
		make_inst(FADDPIPE,FSUB,e,e,t);

		t = get_rotating_register(CMADregs);
		make_inst(FMULPIPE,FMUL,t,b,c);
		make_inst(FADDPIPE,FADD,f,f,t);
	}
}


/**
 * Queues vector complex MADD: (e + if) = (e+if) + (a+ib)*(c+id)
 */
void complex_madd(int e, int a, int c)
{
	if (have_hummer())
	{
		make_inst(SIMDHUMMERPIPE, VMADD_RR_RI, e, a, c, e);
		make_inst(SIMDHUMMERPIPE, VMADD_MII_IR, e, a, c, e);
	}
	else if (have_mic())
	{
		int t1, t2;

		t1 = allocate_tmp_reg(Fregs);
		t2 = allocate_tmp_reg(Fregs);
		make_inst(SIMDMICPIPE, MIC_MUL, t1, a, c);  // RR,II
		make_inst(SIMDMICSWIZPIPE, MIC_MUL_SWIZ, t2, a, c, MIC_SWIZ_CDAB);  // RI,IR
		make_inst(SIMDMICSWIZPIPE, MIC_SUB_MASK_SWIZ, t1, MIC_MASK_ODD,  t1, t1, MIC_SWIZ_CDAB);  // RR-II,??
		make_inst(SIMDMICSWIZPIPE, MIC_ADD_MASK_SWIZ, t1, MIC_MASK_EVEN, t2, t2, MIC_SWIZ_CDAB);  // RR-II,RI+IR
		make_inst(SIMDMICPIPE, MIC_ADD, e, t1, e);
		free_reg(Fregs, t1);
		free_reg(Fregs, t2);
	}
	else
	{
		for (int v = 0; v < nsimd(); v++)
		{
			int f = e+1; /*RE/IM innermost for hummer style SIMD*/
			int b = a+1;
			int d = c+1;
			queue_cmadd(e,f,a,b,c,d);
			e+=2; a+=2; c+=2;
		}
	}
}


/**
 * e += (r*c)
 */
void complex_real_madd(int e, int r, int c)
{
	complex_real_madd(e,r,c,e);
}


/**
 * e = (r*c) + g
 */
void complex_real_madd(int e, int r, int c, int g)
{
	if (have_hummer())
	{
		make_inst(SIMDHUMMERPIPE, VMADD_RR_RI, e, r, c, g);
	}
	else if (have_mic())
	{
		// We assume that r is a vector reg of form (r+ir), i.e., we don't need any broadcast/splat/etc.
		// Must reduce operand count: 4=>3
		if (e == g)      make_inst(SIMDMICPIPE, MIC_FMADD231, e, r, c);
		else if (e == c) make_inst(SIMDMICPIPE, MIC_FMADD213, e, r, g);
		else if (e == r) make_inst(SIMDMICPIPE, MIC_FMADD213, e, c, g);
		else {
			make_inst(SIMDMICPIPE, MIC_MOV, e, g);
			make_inst(SIMDMICPIPE, MIC_FMADD231, e, r, c);
		}

		/*
		// -------------------
		// Assuming we have a vector (r+ix), we must broadcast r somehow to (r+ir)...
		// Unfortunately there's no swizzle for that. But the permute can be used as a copy to solve the 4-operand problem.
		make_inst(SIMDPERMUTEPIPE, VPERMIMM, e, r, MIC_PERM_CCAA);
		make_inst(SIMDMICPIPE, MIC_FMADD213, e, c, g);
		*/
	}
	else
	{
		for (int v = 0; v < nsimd(); v++)
		{
			int f = e+1;
			int h = g+1;
			int d = c+1;

			if (have_fmadd)
			{
				make_inst(FMACPIPE, FMADD, e, r, c, g);
				make_inst(FMACPIPE, FMADD, f, r, d, h);
			}
			else
			{
				int t;
				t = get_rotating_register(CMADregs);
				make_inst(FMULPIPE, FMUL, t, r, c);
				make_inst(FADDPIPE, FADD, e, g, t);
				t = get_rotating_register(CMADregs);
				make_inst(FMULPIPE, FMUL, t, r, d);
				make_inst(FADDPIPE, FADD, f, h, t);
			}
			e+=2; c+=2; g+=2;
		}
	}
}


/**
 * (e+if) = (r+ir) * (c+id)
 */
void complex_real_mul(int e, int r, int c)
{
	if (have_hummer())
	{
		make_inst(SIMDHUMMERPIPE, VMUL_RR_RI, e, r, c, e);
	}
	else if (have_mic())
	{
		// TODO: which implementation is better:
		//  1) mask_mul(r,c) and mask_mul(r{cdab}, c)   lat=4+6
		//  2) shuffle(r) and mul(shuf,c)               lat=6+4

		//make_inst(SIMDPERMUTEPIPE, VPERMIMM, e, r, MIC_PERM_CCAA);
		//make_inst(SIMDMICPIPE, MIC_MUL, e, e, c);

		make_inst(SIMDMICPIPE, MIC_MUL, e, r, c);  // XXX: Can we assume that vector r has the form (r+ir)? Or do we have to do a shuffle?
	}
	else
	{
		for (int v = 0; v < nsimd(); v++)
		{
			int f = e+1;
			int d = c+1;
			make_inst(FMULPIPE,FMUL,e,r,c);
			make_inst(FMULPIPE,FMUL,f,r,d);
			e+=2; c+=2;
		}
	}
}


/**
 * Scalar complex MADD: (e+if) = (e+if) + ((a+ib)*(c+id))
 */
void queue_cmadd(int e,int f,int a,int b,int c,int d)
{
	if (have_fmadd)
	{
		if (have_mic()) {
			make_inst(FMACPIPE, FMADD231, e, a, c);   // e += a*c
			make_inst(FMACPIPE, FMADD231, f, a, d);   // f += a*d
			make_inst(FMACPIPE, FNMADD231, e, b, d);  // e -= b*d
			make_inst(FMACPIPE, FMADD231, f, b, c);   // f += b*c
		}
		else {
			make_inst(FMACPIPE,FMADD,e,a,c,e);
			make_inst(FMACPIPE,FMADD,f,a,d,f);
			make_inst(FMACPIPE,FNMSUB,e,b,d,e);
			make_inst(FMACPIPE,FMADD,f,b,c,f);
		}
	}
	else
	{
		int t;

		t = get_rotating_register(CMADregs);
		make_inst(FMULPIPE,FMUL,t,a,c);
		make_inst(FADDPIPE,FADD,e,e,t);

		t = get_rotating_register(CMADregs);
		make_inst(FMULPIPE,FMUL,t,a,d);
		make_inst(FADDPIPE,FADD,f,f,t);

		t = get_rotating_register(CMADregs);
		make_inst(FMULPIPE,FMUL,t,b,d);
		make_inst(FADDPIPE,FSUB,e,e,t);

		t = get_rotating_register(CMADregs);
		make_inst(FMULPIPE,FMUL,t,b,c);
		make_inst(FADDPIPE,FADD,f,f,t);
	}
}


/**
 * Queues non-destructive vector complex MADD: (h + ig) = (e+if) + (a+ib)*(c+id)
 */
void complex_maddto(int h, int e, int a, int c)
{
	if (have_hummer())
	{
		make_inst(SIMDHUMMERPIPE, VMADD_RR_RI, h, a, c, e);
		make_inst(SIMDHUMMERPIPE, VMADD_MII_IR, h, a, c, h);
	}
	else if (have_mic())
	{
		int t1, t2;

		t1 = allocate_tmp_reg(Fregs);
		t2 = allocate_tmp_reg(Fregs);
		make_inst(SIMDMICPIPE, MIC_MUL, t1, a, c);  // RR,II
		make_inst(SIMDMICSWIZPIPE, MIC_MUL_SWIZ, t2, a, c, MIC_SWIZ_CDAB);  // RI,IR
		make_inst(SIMDMICSWIZPIPE, MIC_SUB_MASK_SWIZ, t1, MIC_MASK_ODD,  t1, t1, MIC_SWIZ_CDAB);  // RR-II,??
		make_inst(SIMDMICSWIZPIPE, MIC_ADD_MASK_SWIZ, t1, MIC_MASK_EVEN, t2, t2, MIC_SWIZ_CDAB);  // RR-II,RI+IR
		make_inst(SIMDMICPIPE, MIC_ADD, h, t1, e);
		free_reg(Fregs, t1);
		free_reg(Fregs, t2);
	}
	else
	{
		for (int v = 0; v < nsimd(); v++)
		{
			int f = e+1;
			int b = a+1;
			int d = c+1;
			int g = h+1;
			queue_cmaddto(h,g,e,f,a,b,c,d);
			e+=2; a+=2; c+=2; h+=2;
		}
	}
}


/**
 * Scalar complex MADD: (h+ig) = (e+if) + ((a+ib)*(c+id))
 */
void queue_cmaddto(int h,int g,int e,int f,int a,int b,int c,int d)
{
	if (have_fmadd) {
		if (have_mic()) {
			make_inst(FMOVPIPE, FMOV, h, e);
			make_inst(FMOVPIPE, FMOV, g, f);
			make_inst(FMACPIPE, FMADD231, h, a, c);   // h += a*c
			make_inst(FMACPIPE, FMADD231, g, a, d);   // g += a*d
			make_inst(FMACPIPE, FNMADD231, h, b, d);  // h -= b*d
			make_inst(FMACPIPE, FMADD231, g, b, c);   // g += b*c
		}
		else {
			make_inst(FMACPIPE,FMADD,h,a,c,e);
			make_inst(FMACPIPE,FMADD,g,a,d,f);
			make_inst(FMACPIPE,FNMSUB,h,b,d,h);
			make_inst(FMACPIPE,FMADD,g,b,c,g);
		}
	} else {
		int t;
/*Want to leave e and f undisturbed so do the copy in the first one*/
/*Otherwise risk different behaviour on FMADD and non-FMADD chips*/
		t = get_rotating_register(CMADregs);
		make_inst(FMULPIPE,FMUL,t,a,c);
		make_inst(FADDPIPE,FADD,h,e,t);

		t = get_rotating_register(CMADregs);
		make_inst(FMULPIPE,FMUL,t,a,d);
		make_inst(FADDPIPE,FADD,h,f,t);

		t = get_rotating_register(CMADregs);
		make_inst(FMULPIPE,FMUL,t,b,d);
		make_inst(FADDPIPE,FSUB,h,h,t);

		t = get_rotating_register(CMADregs);
		make_inst(FMULPIPE,FMUL,t,b,c);
		make_inst(FADDPIPE,FADD,g,g,t);
	}
}


/**
 * Queues vector complex:
 * (e + if) = conj(a+ib) * (c+id)
 *          = (a-ib) * (c+id)
 *          = (ac+bd) + i(ad-bc)
 */
void complex_conjmul(int e, int a, int c)
{
	if (have_hummer())
	{
		// complex_constants_assert(cmplx_zero);
		make_inst(SIMDHUMMERPIPE, VMUL_RR_RI, e, a, c, reg_tmp);
		make_inst(SIMDHUMMERPIPE, VMADD_II_MIR, e, a, c, e);
	}
	else if (have_mic())
	{
		int t;

		t = allocate_tmp_reg(Fregs);
		make_inst(SIMDMICPIPE, MIC_MUL, e, a, c);  // RR,II
		make_inst(SIMDMICSWIZPIPE, MIC_MUL_SWIZ, t, a, c, MIC_SWIZ_CDAB);  // RI,IR
		make_inst(SIMDMICSWIZPIPE, MIC_ADD_MASK_SWIZ, e, MIC_MASK_ODD,  e, e, MIC_SWIZ_CDAB);  // RR+II,??
		make_inst(SIMDMICSWIZPIPE, MIC_SUBR_MASK_SWIZ, e, MIC_MASK_EVEN, t, t, MIC_SWIZ_CDAB);  // RR+II,RI-IR
		free_reg(Fregs, t);
	}
	else
	{
		for (int v=0; v<nsimd(); v++)
		{
			int f = e+1;
			int b = a+1;
			int d = c+1;
			queue_conjmul(e,f,a,b,c,d);
			e+=2; a+=2; c+=2;
		}
	}
}


/**
 * (e+if) = conj((a+ib))*(c+id)
 *        = (ac+bd)+i(ad-bc)
 */
void queue_conjmul(int e, int f, int a, int b, int c, int d)
{
	if (have_fmadd)
	{
		if (have_mic()) {
			make_inst(FMULPIPE, FMUL, e, a, c);
			make_inst(FMULPIPE, FMUL, f, a, d);
			make_inst(FMACPIPE, FMADD231, e, b, d);   // e = ac+bd
			make_inst(FMACPIPE, FNMADD231, f, b, c);  // f = ad-bc
		}
		else {
			make_inst(FMULPIPE,FMUL, e, a, c);        // e = a*c
			make_inst(FMULPIPE,FMUL, f, a, d);        // f = a*d
			make_inst(FMACPIPE,FMADD, e, b, d, e);    // e = (b*d)+e = ac+bd
			make_inst(FMACPIPE,FNMSUB, f, b, c, f);   // f = f-(b*c) = ad-bc
		}
	}
	else
	{
		// XXX: QUESTION: is there a bug? "if" and "else" block seem to do different things.....?
		int t;

		make_inst(FMULPIPE,FMUL,e,a,c);           // e = a*c
		make_inst(FMULPIPE,FMUL,f,a,d);           // f = a*d

		t = get_rotating_register(CMADregs);
		make_inst(FMULPIPE,FMUL,t,b,c);          // t = b*c
		make_inst(FADDPIPE,FADD,e,e,t);          // e = e+t = ac+bc

		t = get_rotating_register(CMADregs);
		make_inst(FMULPIPE,FMUL,t,b,d);          // t = b*d
		make_inst(FADDPIPE,FSUB,f,f,t);          // f = f-t = ad-bd
	}
}


/**
 * Vector complex:
 * e = e + (conj(a) * c)
 *   = (e+if) + (ad+bc) + i(ac-bd)
 */
void complex_conjmadd(int e,int a,int c)
{
	if (have_hummer())
	{
		make_inst(SIMDHUMMERPIPE,VMADD_RR_RI,e,a,c,e);
		make_inst(SIMDHUMMERPIPE,VMADD_II_MIR,e,a,c,e);
	}
	else if (have_mic())
	{
		int t1, t2;

		t1 = allocate_tmp_reg(Fregs);
		t2 = allocate_tmp_reg(Fregs);
		make_inst(SIMDMICPIPE, MIC_MUL, t1, a, c);  // RR,II
		make_inst(SIMDMICSWIZPIPE, MIC_MUL_SWIZ, t2, a, c, MIC_SWIZ_CDAB);  // RI,IR
		make_inst(SIMDMICSWIZPIPE, MIC_ADD_MASK_SWIZ, t1, MIC_MASK_ODD,  t1, t1, MIC_SWIZ_CDAB);  // RR+II,??
		make_inst(SIMDMICSWIZPIPE, MIC_SUBR_MASK_SWIZ, t1, MIC_MASK_EVEN, t2, t2, MIC_SWIZ_CDAB);  // RR+II,RI-IR
		make_inst(SIMDMICPIPE, MIC_ADD, e, t1, e);
		free_reg(Fregs, t1);
		free_reg(Fregs, t2);
	}
	else
	{
		for (int v=0; v<nsimd(); v++)
		{
			int f = e+1;
			int b = a+1;
			int d = c+1;
			queue_conjmadd(e,f,a,b,c,d);
			e+=2; a+=2; c+=2;
		}
	}
}

/**
 * (e+if) = (e+if) + (conj(a+ib) * (c+id))
 *        = (e+if) + ((a-ib) * (c+id))
 *        = (e+if) + (ac+bd) + i(ad-bc)
 */
void queue_conjmadd(int e,int f,int a,int b,int c,int d)
{
	if (have_fmadd) {
		if (have_mic()) {
			make_inst(FMACPIPE, FMADD231, e, a, c);
			make_inst(FMACPIPE, FMADD231, f, a, d);
			make_inst(FMACPIPE, FMADD231, e, b, d);   // e = e_org + ac + bd
			make_inst(FMACPIPE, FNMADD231, f, b, c);  // f = f_org + ad - bc
		}
		else {
			make_inst(FMACPIPE,FMADD,e,a,c,e);   // e += ac
			make_inst(FMACPIPE,FMADD,f,a,d,f);   // f += ad
			make_inst(FMACPIPE,FMADD,e,b,d,e);   // e = e_org + ac + bd
			make_inst(FMACPIPE,FNMSUB,f,b,c,f);  // f = f_org + ad - bc
		}
	}
	else {
		int t;
		t = get_rotating_register(CMADregs);
		make_inst(FMULPIPE,FMUL,t,a,c);
		make_inst(FADDPIPE,FADD,e,e,t);

		t = get_rotating_register(CMADregs);
		make_inst(FMULPIPE,FMUL,t,a,d);
		make_inst(FADDPIPE,FADD,f,f,t);

		t = get_rotating_register(CMADregs);
		make_inst(FMULPIPE,FMUL,t,b,d);
		make_inst(FADDPIPE,FADD,e,e,t);

		t = get_rotating_register(CMADregs);
		make_inst(FMULPIPE,FMUL,t,b,c);
		make_inst(FADDPIPE,FSUB,f,f,t);
	}
}

/**
 * h = e + conj(a)*c
 */
void complex_conjmaddto(int h, int e,int a,int c)
{
	if (have_hummer())
	{
		make_inst(SIMDHUMMERPIPE,VMADD_RR_RI,h,a,c,e);
		make_inst(SIMDHUMMERPIPE,VMADD_II_MIR,h,a,c,h);
	}
	else if (have_mic())
	{
		int t1, t2;

		t1 = allocate_tmp_reg(Fregs);
		t2 = allocate_tmp_reg(Fregs);
		make_inst(SIMDMICPIPE, MIC_MUL, t1, a, c);  // RR,II
		make_inst(SIMDMICSWIZPIPE, MIC_MUL_SWIZ, t2, a, c, MIC_SWIZ_CDAB);  // RI,IR
		make_inst(SIMDMICSWIZPIPE, MIC_ADD_MASK_SWIZ, t1, MIC_MASK_ODD,  t1, t1, MIC_SWIZ_CDAB);  // RR+II,??
		make_inst(SIMDMICSWIZPIPE, MIC_SUBR_MASK_SWIZ, t1, MIC_MASK_EVEN, t2, t2, MIC_SWIZ_CDAB);  // RR+II,RI-IR
		make_inst(SIMDMICPIPE, MIC_ADD, h, t1, e);
		free_reg(Fregs, t1);
		free_reg(Fregs, t2);
	}
	else
	{
		for (int v = 0; v<nsimd(); v++)
		{
			int f = e+1;
			int b = a+1;
			int d = c+1;
			int g = h+1;
			queue_conjmaddto(h,g,e,f,a,b,c,d);
			e+=2; a+=2; c+=2; h+=2;
		}
	}
}

void queue_conjmaddto(int h,int g,int e,int f,int a,int b,int c,int d)
{
	if (have_fmadd) {
		if (have_mic()) {
			make_inst(FMOVPIPE, FMOV, h, e);
			make_inst(FMOVPIPE, FMOV, g, f);
			make_inst(FMACPIPE, FMADD231, h, a, c);
			make_inst(FMACPIPE, FMADD231, g, a, d);
			make_inst(FMACPIPE, FMADD231, h, b, d);   // h = e + ac + bd
			make_inst(FMACPIPE, FNMADD231, g, b, c);  // g = f + ad - bc
		}
		else {
			make_inst(FMACPIPE,FMADD,h,a,c,e);   // h = ac + e
			make_inst(FMACPIPE,FMADD,g,a,d,f);   // g = ad + f
			make_inst(FMACPIPE,FMADD,h,b,d,h);   // h = h + bd = e + ac + bd
			make_inst(FMACPIPE,FNMSUB,g,b,c,g);  // g = g - bc = f + ad - bc
		}
	}
	else {
		int t;
		t = get_rotating_register(CMADregs);
		make_inst(FMULPIPE,FMUL,t,a,c);
		make_inst(FADDPIPE,FADD,h,e,t);

		t = get_rotating_register(CMADregs);
		make_inst(FMULPIPE,FMUL,t,a,d);
		make_inst(FADDPIPE,FADD,g,f,t);

		t = get_rotating_register(CMADregs);
		make_inst(FMULPIPE,FMUL,t,b,d);
		make_inst(FADDPIPE,FADD,h,h,t);

		t = get_rotating_register(CMADregs);
		make_inst(FMULPIPE,FMUL,t,b,c);
		make_inst(FADDPIPE,FSUB,g,g,t);
	}
}


void complex_six_cmuls(
	int e1,int a1,int c1,
	int e2,int a2,int c2,
	int e3,int a3,int c3,
	int e4,int a4,int c4,
	int e5,int a5,int c5,
	int e6,int a6,int c6)
{
	if (have_hummer())
	{
		//complex_constants_assert(cmplx_zero);
		make_inst(SIMDHUMMERPIPE,VMUL_RR_RI,e1,a1,c1,reg_tmp);
		make_inst(SIMDHUMMERPIPE,VMUL_RR_RI,e2,a2,c2,reg_tmp);
		make_inst(SIMDHUMMERPIPE,VMUL_RR_RI,e3,a3,c3,reg_tmp);
		make_inst(SIMDHUMMERPIPE,VMUL_RR_RI,e4,a4,c4,reg_tmp);
		make_inst(SIMDHUMMERPIPE,VMUL_RR_RI,e5,a5,c5,reg_tmp);
		make_inst(SIMDHUMMERPIPE,VMUL_RR_RI,e6,a6,c6,reg_tmp);

		make_inst(SIMDHUMMERPIPE,VMADD_MII_IR,e1,a1,c1,e1);
		make_inst(SIMDHUMMERPIPE,VMADD_MII_IR,e2,a2,c2,e2);
		make_inst(SIMDHUMMERPIPE,VMADD_MII_IR,e3,a3,c3,e3);
		make_inst(SIMDHUMMERPIPE,VMADD_MII_IR,e4,a4,c4,e4);
		make_inst(SIMDHUMMERPIPE,VMADD_MII_IR,e5,a5,c5,e5);
		make_inst(SIMDHUMMERPIPE,VMADD_MII_IR,e6,a6,c6,e6);
	}
	else if (have_mic())
	{
		// TODO: maybe fine-tune by hand?
		complex_mul(e1, a1, c1);
		complex_mul(e2, a2, c2);
		complex_mul(e3, a3, c3);
		complex_mul(e4, a4, c4);
		complex_mul(e5, a5, c5);
		complex_mul(e6, a6, c6);
	}
	else
	{
		complex_three_cmuls(e1,a1,c1,
			e2,a2,c2,
			e3,a3,c3);
		complex_three_cmuls(e4,a4,c4,
			e5,a5,c5,
			e6,a6,c6);
	}
}
void complex_twospin_color_dot(int u1,int u2,int u3,
			       int in_h1,int in_h2,int in_h3,
			       int in_l1,int in_l2,int in_l3,
			       int out_h,int out_l)
{
  if ( have_mic() ) {
    int t_l, t_h;
    t_h = allocate_tmp_reg(Fregs);
    t_l = allocate_tmp_reg(Fregs);
    make_inst(SIMDMICPIPE, MIC_MUL, out_h, u1, in_h1);  // RR,II
    make_inst(SIMDMICSWIZPIPE, MIC_MUL_SWIZ, t_h, u1, in_h1,MIC_SWIZ_CDAB); //RI,IR
    make_inst(SIMDMICPIPE, MIC_MUL, out_l, u1, in_l1);  // RR,II
    make_inst(SIMDMICSWIZPIPE, MIC_MUL_SWIZ, t_l, u1, in_l1,MIC_SWIZ_CDAB); //RI,IR

    make_inst(SIMDMICPIPE, MIC_FMADD231, out_h, u2, in_h2);  // RR,II
    make_inst(SIMDMICSWIZPIPE, MIC_FMADD231_SWIZ, t_h, u2, in_h2,MIC_SWIZ_CDAB); //RI,IR
    make_inst(SIMDMICPIPE, MIC_FMADD231, out_l, u2, in_l2);  // RR,II
    make_inst(SIMDMICSWIZPIPE, MIC_FMADD231_SWIZ, t_l, u2, in_l2,MIC_SWIZ_CDAB); //RI,IR

    make_inst(SIMDMICPIPE, MIC_FMADD231, out_h, u3, in_h3);  // RR,II
    make_inst(SIMDMICSWIZPIPE, MIC_FMADD231_SWIZ, t_h, u3, in_h3,MIC_SWIZ_CDAB); //RI,IR
    make_inst(SIMDMICPIPE, MIC_FMADD231, out_l, u3, in_l3);  // RR,II
    make_inst(SIMDMICSWIZPIPE, MIC_FMADD231_SWIZ, t_l, u3, in_l3,MIC_SWIZ_CDAB); //RI,IR

    make_inst(SIMDMICSWIZPIPE, MIC_SUB_MASK_SWIZ, out_h, MIC_MASK_ODD, out_h, out_h, MIC_SWIZ_CDAB);  // RR-II,??
    make_inst(SIMDMICSWIZPIPE, MIC_ADD_MASK_SWIZ, out_h, MIC_MASK_EVEN, t_h, t_h, MIC_SWIZ_CDAB);  // RR+II,RI-IR
    make_inst(SIMDMICSWIZPIPE, MIC_SUB_MASK_SWIZ, out_l, MIC_MASK_ODD, out_l, out_l, MIC_SWIZ_CDAB);  // RR-II,??
    make_inst(SIMDMICSWIZPIPE, MIC_ADD_MASK_SWIZ, out_l, MIC_MASK_EVEN, t_l, t_l, MIC_SWIZ_CDAB);  // RR+II,RI-IR
    free_reg(Fregs, t_l);
    free_reg(Fregs, t_h);
  } else {
    complex_mul(out_h,u1,in_h1);
    complex_mul(out_l,u1,in_l1);
    complex_madd(out_h,u2,in_h2);
    complex_madd(out_l,u2,in_l2);
    complex_madd(out_h,u3,in_h3);
    complex_madd(out_l,u3,in_l3);
  }
}

void complex_six_cmaddtos(
		       int h1, int e1,int a1,int c1,
                       int h2, int e2,int a2,int c2,
                       int h3, int e3,int a3,int c3,
		       int h4, int e4,int a4,int c4,
                       int h5, int e5,int a5,int c5,
                       int h6, int e6,int a6,int c6)
{
  if ( have_hummer() ){

    //complex_constants_assert(cmplx_zero);
    make_inst(SIMDHUMMERPIPE,VMADD_RR_RI,e1,a1,c1,e1);
    make_inst(SIMDHUMMERPIPE,VMADD_RR_RI,e2,a2,c2,e2);
    make_inst(SIMDHUMMERPIPE,VMADD_RR_RI,e3,a3,c3,e3);
    make_inst(SIMDHUMMERPIPE,VMADD_RR_RI,e4,a4,c4,e4);
    make_inst(SIMDHUMMERPIPE,VMADD_RR_RI,e5,a5,c5,e5);
    make_inst(SIMDHUMMERPIPE,VMADD_RR_RI,e6,a6,c6,e6);

    make_inst(SIMDHUMMERPIPE,VMADD_MII_IR,h1,a1,c1,e1);// t1 & t5
    make_inst(SIMDHUMMERPIPE,VMADD_MII_IR,h2,a2,c2,e2);// t2 & t6
    make_inst(SIMDHUMMERPIPE,VMADD_MII_IR,h3,a3,c3,e3);// t3 & t0
    make_inst(SIMDHUMMERPIPE,VMADD_MII_IR,h4,a4,c4,e4);// t4 & t5
    make_inst(SIMDHUMMERPIPE,VMADD_MII_IR,h5,a5,c5,e5);// t5 & t6
    make_inst(SIMDHUMMERPIPE,VMADD_MII_IR,h6,a6,c6,e6);// t6 & t0

  } else { 
    /*
    complex_three_cmaddtos(
			   h1, e1, a1, c1,
			   h2, e2, a2, c2,
			   h3, e3, a3, c3);
    complex_three_cmaddtos(
			   h4, e4, a4, c4,
			   h5, e5, a5, c5,
			   h6, e6, a6, c6);
    */
    complex_maddto(h1, e1, a1, c1);
    complex_maddto(h2, e2, a2, c2);
    complex_maddto(h3, e3, a3, c3);
    complex_maddto(h4, e4, a4, c4);
    complex_maddto(h5, e5, a5, c5);
    complex_maddto(h6, e6, a6, c6);
  }
}


void complex_six_cmadds(int e1,int a1,int c1,
	int e2,int a2,int c2,
	int e3,int a3,int c3,
	int e4,int a4,int c4,
	int e5,int a5,int c5,
	int e6,int a6,int c6)
{
	if (have_hummer())
	{
		//complex_constants_assert(cmplx_zero);
		make_inst(SIMDHUMMERPIPE,VMADD_RR_RI,e1,a1,c1,e1);
		make_inst(SIMDHUMMERPIPE,VMADD_RR_RI,e2,a2,c2,e2);
		make_inst(SIMDHUMMERPIPE,VMADD_RR_RI,e3,a3,c3,e3);
		make_inst(SIMDHUMMERPIPE,VMADD_RR_RI,e4,a4,c4,e4);
		make_inst(SIMDHUMMERPIPE,VMADD_RR_RI,e5,a5,c5,e5);
		make_inst(SIMDHUMMERPIPE,VMADD_RR_RI,e6,a6,c6,e6);

		make_inst(SIMDHUMMERPIPE,VMADD_MII_IR,e1,a1,c1,e1);
		make_inst(SIMDHUMMERPIPE,VMADD_MII_IR,e2,a2,c2,e2);
		make_inst(SIMDHUMMERPIPE,VMADD_MII_IR,e3,a3,c3,e3);
		make_inst(SIMDHUMMERPIPE,VMADD_MII_IR,e4,a4,c4,e4);
		make_inst(SIMDHUMMERPIPE,VMADD_MII_IR,e5,a5,c5,e5);
		make_inst(SIMDHUMMERPIPE,VMADD_MII_IR,e6,a6,c6,e6);
	}
	else
	{
	  complex_three_cmadds(e1,a1,c1,
			       e2,a2,c2,
			       e3,a3,c3);
	  complex_three_cmadds(e4,a4,c4,
			       e5,a5,c5,
			       e6,a6,c6);
	}
}




void complex_three_cmuls(
	int e1,int a1,int c1,
	int e2,int a2,int c2,
	int e3,int a3,int c3)
{
	if (have_hummer())
	{
		//complex_constants_assert(cmplx_zero);
		make_inst(SIMDHUMMERPIPE,VMUL_RR_RI,e1,a1,c1,reg_tmp);
		make_inst(SIMDHUMMERPIPE,VMUL_RR_RI,e2,a2,c2,reg_tmp);
		make_inst(SIMDHUMMERPIPE,VMUL_RR_RI,e3,a3,c3,reg_tmp);

		make_inst(SIMDHUMMERPIPE,VMADD_MII_IR,e1,a1,c1,e1);
		make_inst(SIMDHUMMERPIPE,VMADD_MII_IR,e2,a2,c2,e2);
		make_inst(SIMDHUMMERPIPE,VMADD_MII_IR,e3,a3,c3,e3);
	/*
	* } else if (have_fmadd){
	*
	* int f1 = e1 + 1;
	* int f2 = e2 + 1;
	* int f3 = e3 + 1;

	* int b1 = a1 + 1;
	* int b2 = a2 + 1;
	* int b3 = a3 + 1;

	* int d1 = c1 + 1;
	* int d2 = c2 + 1;
	*  int d3 = c3 + 1;
	*    queue_three_cmuls(e1, f1, a1, b1, c1, d1,
	*		      e2, f2, a2, b2, c2, d2,
	*		      e3, f3, a3, b3, c3, d3);
	*/
	}
	else if (have_mic())
	{
		int t1, t2, t3;

		t1 = allocate_tmp_reg(Fregs);
		t2 = allocate_tmp_reg(Fregs);
		t3 = allocate_tmp_reg(Fregs);

		make_inst(SIMDMICPIPE, MIC_MUL, e1, a1, c1);
		make_inst(SIMDMICPIPE, MIC_MUL, e2, a2, c2);
		make_inst(SIMDMICPIPE, MIC_MUL, e3, a3, c3);

		make_inst(SIMDMICSWIZPIPE, MIC_MUL_SWIZ, t1, a1, c1, MIC_SWIZ_CDAB);
		make_inst(SIMDMICSWIZPIPE, MIC_MUL_SWIZ, t2, a2, c2, MIC_SWIZ_CDAB);
		make_inst(SIMDMICSWIZPIPE, MIC_MUL_SWIZ, t3, a3, c3, MIC_SWIZ_CDAB);

		make_inst(SIMDMICSWIZPIPE, MIC_SUB_MASK_SWIZ, e1, MIC_MASK_ODD,  e1, e1, MIC_SWIZ_CDAB);
		make_inst(SIMDMICSWIZPIPE, MIC_SUB_MASK_SWIZ, e2, MIC_MASK_ODD,  e2, e2, MIC_SWIZ_CDAB);
		make_inst(SIMDMICSWIZPIPE, MIC_SUB_MASK_SWIZ, e3, MIC_MASK_ODD,  e3, e3, MIC_SWIZ_CDAB);

		make_inst(SIMDMICSWIZPIPE, MIC_ADD_MASK_SWIZ, e1, MIC_MASK_EVEN, t1, t1, MIC_SWIZ_CDAB);
		make_inst(SIMDMICSWIZPIPE, MIC_ADD_MASK_SWIZ, e2, MIC_MASK_EVEN, t2, t2, MIC_SWIZ_CDAB);
		make_inst(SIMDMICSWIZPIPE, MIC_ADD_MASK_SWIZ, e3, MIC_MASK_EVEN, t3, t3, MIC_SWIZ_CDAB);

		free_reg(Fregs, t1);
		free_reg(Fregs, t2);
		free_reg(Fregs, t3);
	}
	else
	{
		complex_mul(e1, a1, c1);
		complex_mul(e2, a2, c2);
		complex_mul(e3, a3, c3);
	}
}


void queue_three_cmuls(
	int e1, int f1, int a1, int b1, int c1, int d1,
	int e2, int f2, int a2, int b2, int c2, int d2,
	int e3, int f3, int a3, int b3, int c3, int d3)
{
	queue_cmul(e1,f1,a1,b1,c1,d1);
	queue_cmul(e2,f2,a2,b2,c2,d2);
	queue_cmul(e3,f3,a3,b3,c3,d3);
}


/**
 * "Destructive" complex multiply-add:
 * e1 = (a1 * c1) + e1
 * e2 = ...
 * e3 = ...
 */
void complex_three_cmadds(
	int e1,int a1,int c1,
	int e2,int a2,int c2,
	int e3,int a3,int c3)
{
	if (have_hummer())
	{
		make_inst(SIMDHUMMERPIPE,VMADD_RR_RI,e1,a1,c1,e1);
		make_inst(SIMDHUMMERPIPE,VMADD_RR_RI,e2,a2,c2,e2);
		make_inst(SIMDHUMMERPIPE,VMADD_RR_RI,e3,a3,c3,e3);

		make_inst(SIMDHUMMERPIPE,VMADD_MII_IR,e1,a1,c1,e1);
		make_inst(SIMDHUMMERPIPE,VMADD_MII_IR,e2,a2,c2,e2);
		make_inst(SIMDHUMMERPIPE,VMADD_MII_IR,e3,a3,c3,e3);
	}
	else if (have_mic() && 0)
	{
		// TODO: check if this works correctly!
		int t1, t2, t3, t4, t5, t6;

		t1 = allocate_tmp_reg(Fregs);
		t2 = allocate_tmp_reg(Fregs);
		t3 = allocate_tmp_reg(Fregs);
		t4 = allocate_tmp_reg(Fregs);
		t5 = allocate_tmp_reg(Fregs);

		// kinda turbulent code. Want to avoid stalls and keep no. of temp. regs small.
		make_inst(SIMDMICPIPE, MIC_MUL, t1, a1, c1);
		make_inst(SIMDMICPIPE, MIC_MUL, t2, a2, c2);
		make_inst(SIMDMICPIPE, MIC_MUL, t3, a3, c3);

		make_inst(SIMDMICSWIZPIPE, MIC_SUB_MASK_SWIZ, t1, MIC_MASK_ODD,  t1, t1, MIC_SWIZ_CDAB);
		make_inst(SIMDMICSWIZPIPE, MIC_MUL_SWIZ, t4, a1, c1, MIC_SWIZ_CDAB);
		make_inst(SIMDMICSWIZPIPE, MIC_SUB_MASK_SWIZ, t2, MIC_MASK_ODD,  t2, t2, MIC_SWIZ_CDAB);
		make_inst(SIMDMICSWIZPIPE, MIC_MUL_SWIZ, t5, a2, c2, MIC_SWIZ_CDAB);
		make_inst(SIMDMICSWIZPIPE, MIC_SUB_MASK_SWIZ, t3, MIC_MASK_ODD,  t3, t3, MIC_SWIZ_CDAB);

		make_inst(SIMDMICSWIZPIPE, MIC_ADD_MASK_SWIZ, t1, MIC_MASK_EVEN, t4, t4, MIC_SWIZ_CDAB);  // t4 now free

		t6 = t4;  // re-use t4
		make_inst(SIMDMICSWIZPIPE, MIC_MUL_SWIZ, t6, a3, c3, MIC_SWIZ_CDAB);

		make_inst(SIMDMICSWIZPIPE, MIC_ADD_MASK_SWIZ, t2, MIC_MASK_EVEN, t5, t5, MIC_SWIZ_CDAB);
		make_inst(SIMDMICSWIZPIPE, MIC_ADD_MASK_SWIZ, t3, MIC_MASK_EVEN, t6, t6, MIC_SWIZ_CDAB);

		make_inst(SIMDMICPIPE, MIC_ADD, e1, t1, e1);
		make_inst(SIMDMICPIPE, MIC_ADD, e2, t2, e2);
		make_inst(SIMDMICPIPE, MIC_ADD, e3, t3, e3);

		free_reg(Fregs, t1);
		free_reg(Fregs, t2);
		free_reg(Fregs, t3);
		free_reg(Fregs, t4);
		free_reg(Fregs, t5);
	}
	else
	{
		complex_madd(e1,a1,c1);
		complex_madd(e2,a2,c2);
		complex_madd(e3,a3,c3);
	}
}


void queue_three_cmadds(int e1, int f1, int a1, int b1, int c1, int d1,
	int e2, int f2, int a2, int b2, int c2, int d2,
	int e3, int f3, int a3, int b3, int c3, int d3)
{
	if (have_fmadd)
	{
		if (have_mic()) {
			// QUESTION: might be ok to simply call queue_cmadd(), which has a MIC-specific implementation, and let the scheduler re-order instructions?
			make_inst(FMACPIPE, FMADD231, e1, a1, c1);
			make_inst(FMACPIPE, FMADD231, f1, a1, d1);
			make_inst(FMACPIPE, FMADD231, e2, a2, c2);
			make_inst(FMACPIPE, FMADD231, f2, a2, d2);
			make_inst(FMACPIPE, FMADD231, e3, a3, c3);
			make_inst(FMACPIPE, FMADD231, f3, a3, d3);

			make_inst(FMACPIPE, FNMADD231, e1, b1, d1);
			make_inst(FMACPIPE, FMADD231,  f1, b1, c1);
			make_inst(FMACPIPE, FNMADD231, e2, b2, d2);
			make_inst(FMACPIPE, FMADD231,  f2, b2, c2);
			make_inst(FMACPIPE, FNMADD231, e3, b3, d3);
			make_inst(FMACPIPE, FMADD231,  f3, b3, c3);
		}
		else {
			make_inst(FMACPIPE,FMADD,e1,a1,c1,e1);
			make_inst(FMACPIPE,FMADD,f1,a1,d1,f1);

			make_inst(FMACPIPE,FMADD,e2,a2,c2,e2);
			make_inst(FMACPIPE,FMADD,f2,a2,d2,f2);

			make_inst(FMACPIPE,FMADD,e3,a3,c3,e3);
			make_inst(FMACPIPE,FMADD,f3,a3,d3,f3);

			make_inst(FMACPIPE,FNMSUB,e1,b1,d1,e1);
			make_inst(FMACPIPE,FMADD,f1,b1,c1,f1);

			make_inst(FMACPIPE,FNMSUB,e2,b2,d2,e2);
			make_inst(FMACPIPE,FMADD,f2,b2,c2,f2);

			make_inst(FMACPIPE,FNMSUB,e3,b3,d3,e3);
			make_inst(FMACPIPE,FMADD,f3,b3,c3,f3);
		}
	}
	else
	{
		queue_cmadd(e1,f1,a1,b1,c1,d1);
		queue_cmadd(e2,f2,a2,b2,c2,d2);
		queue_cmadd(e3,f3,a3,b3,c3,d3);
	}
}

void complex_three_cmaddtos(
	int h1,int e1,int a1,int c1,
	int h2,int e2,int a2,int c2,
	int h3,int e3,int a3,int c3)
{
	if (have_hummer())
	{
		make_inst(SIMDHUMMERPIPE,VMADD_RR_RI,h1,a1,c1,e1);
		make_inst(SIMDHUMMERPIPE,VMADD_RR_RI,h2,a2,c2,e2);
		make_inst(SIMDHUMMERPIPE,VMADD_RR_RI,h3,a3,c3,e3);

		make_inst(SIMDHUMMERPIPE,VMADD_MII_IR,h1,a1,c1,h1);
		make_inst(SIMDHUMMERPIPE,VMADD_MII_IR,h2,a2,c2,h2);
		make_inst(SIMDHUMMERPIPE,VMADD_MII_IR,h3,a3,c3,h3);
	}
	else if (have_mic())
	{
		// TODO: maybe fine-tune
		complex_maddto(h1, e1, a1, c1);
		complex_maddto(h2, e2, a2, c2);
		complex_maddto(h3, e3, a3, c3);
	}
	else
	{
		complex_maddto(h1,e1,a1,c1);
		complex_maddto(h2,e2,a2,c2);
		complex_maddto(h3,e3,a3,c3);
	}
}


void queue_three_cmaddtos(int h1,int g1,int e1,int f1,int a1,int b1,int c1,int d1,
	int h2,int g2,int e2,int f2,int a2,int b2,int c2,int d2,
	int h3,int g3,int e3,int f3,int a3,int b3,int c3,int d3)
{
	if (have_fmadd)
	{
		if (have_mic()) {
			make_inst(FMOVPIPE, FMOV, h1, e1);
			make_inst(FMOVPIPE, FMOV, g1, f1);
			make_inst(FMOVPIPE, FMOV, h2, e2);
			make_inst(FMOVPIPE, FMOV, g2, f2);
			make_inst(FMOVPIPE, FMOV, h3, e3);
			make_inst(FMOVPIPE, FMOV, g3, f3);

			make_inst(FMACPIPE, FMADD231, h1, a1, c1);
			make_inst(FMACPIPE, FMADD231, g1, a1, d1);
			make_inst(FMACPIPE, FMADD231, h2, a2, c2);
			make_inst(FMACPIPE, FMADD231, g2, a2, d2);
			make_inst(FMACPIPE, FMADD231, h3, a3, c3);
			make_inst(FMACPIPE, FMADD231, g3, a3, d3);

			make_inst(FMACPIPE, FNMADD231, h1, b1, d1);
			make_inst(FMACPIPE, FMADD231,  g1, b1, c1);
			make_inst(FMACPIPE, FNMADD231, h2, b2, d2);
			make_inst(FMACPIPE, FMADD231,  g2, b2, c2);
			make_inst(FMACPIPE, FNMADD231, h3, b3, d3);
			make_inst(FMACPIPE, FMADD231,  g3, b3, c3);
		}
		else {
			make_inst(FMACPIPE,FMADD,h1,a1,c1,e1);
			make_inst(FMACPIPE,FMADD,g1,a1,d1,f1);

			make_inst(FMACPIPE,FMADD,h2,a2,c2,e2);
			make_inst(FMACPIPE,FMADD,g2,a2,d2,f2);

			make_inst(FMACPIPE,FMADD,h3,a3,c3,e3);
			make_inst(FMACPIPE,FMADD,g3,a3,d3,f3);

			make_inst(FMACPIPE,FNMSUB,h1,b1,d1,h1);
			make_inst(FMACPIPE,FMADD,g1,b1,c1,g1);

			make_inst(FMACPIPE,FNMSUB,h2,b2,d2,h2);
			make_inst(FMACPIPE,FMADD,g2,b2,c2,g2);

			make_inst(FMACPIPE,FNMSUB,h3,b3,d3,h3);
			make_inst(FMACPIPE,FMADD,g3,b3,c3,g3);
		}
	} else {
		queue_cmaddto(h1,g1,e1,f1,a1,b1,c1,d1);
		queue_cmaddto(h2,g2,e2,f2,a2,b2,c2,d2);
		queue_cmaddto(h3,g3,e3,f3,a3,b3,c3,d3);
	}
}

void complex_three_conjmadds(int e1,int a1,int c1,
	int e2,int a2,int c2,
	int e3,int a3,int c3)
{
	if (have_hummer())
	{
		make_inst(SIMDHUMMERPIPE,VMADD_RR_RI,e1,a1,c1,e1);
		make_inst(SIMDHUMMERPIPE,VMADD_RR_RI,e2,a2,c2,e2);
		make_inst(SIMDHUMMERPIPE,VMADD_RR_RI,e3,a3,c3,e3);

		make_inst(SIMDHUMMERPIPE,VMADD_II_MIR,e1,a1,c1,e1);
		make_inst(SIMDHUMMERPIPE,VMADD_II_MIR,e2,a2,c2,e2);
		make_inst(SIMDHUMMERPIPE,VMADD_II_MIR,e3,a3,c3,e3);
	}
	else if (have_mic())
	{
		// TODO
		complex_conjmadd(e1, a1, c1);
		complex_conjmadd(e2, a2, c2);
		complex_conjmadd(e3, a3, c3);
	}
	else
	{
		complex_conjmadd(e1,a1,c1);
		complex_conjmadd(e2,a2,c2);
		complex_conjmadd(e3,a3,c3);
	}
}

void queue_three_conjmadds(int e1,int f1,int a1,int b1,int c1,int d1,
	int e2,int f2,int a2,int b2,int c2,int d2,
	int e3,int f3,int a3,int b3,int c3,int d3)
{
	if (have_fmadd) {
		if (have_mic()) {
			make_inst(FMACPIPE, FMADD231, e1, a1, c1);
			make_inst(FMACPIPE, FMADD231, f1, a1, d1);
			make_inst(FMACPIPE, FMADD231, e2, a2, c2);
			make_inst(FMACPIPE, FMADD231, f2, a2, d2);
			make_inst(FMACPIPE, FMADD231, e3, a3, c3);
			make_inst(FMACPIPE, FMADD231, f3, a3, d3);

			make_inst(FMACPIPE, FMADD231,  e1, b1, d1);
			make_inst(FMACPIPE, FNMADD231, f1, b1, c1);
			make_inst(FMACPIPE, FMADD231,  e2, b2, d2);
			make_inst(FMACPIPE, FNMADD231, f2, b2, c2);
			make_inst(FMACPIPE, FMADD231,  e3, b3, d3);
			make_inst(FMACPIPE, FNMADD231, f3, b3, c3);
		}
		else {
			make_inst(FMACPIPE,FMADD,e1,a1,c1,e1);
			make_inst(FMACPIPE,FMADD,f1,a1,d1,f1);

			make_inst(FMACPIPE,FMADD,e2,a2,c2,e2);
			make_inst(FMACPIPE,FMADD,f2,a2,d2,f2);

			make_inst(FMACPIPE,FMADD,e3,a3,c3,e3);
			make_inst(FMACPIPE,FMADD,f3,a3,d3,f3);

			make_inst(FMACPIPE,FMADD,e1,b1,d1,e1);
			make_inst(FMACPIPE,FNMSUB,f1,b1,c1,f1);

			make_inst(FMACPIPE,FMADD,e2,b2,d2,e2);
			make_inst(FMACPIPE,FNMSUB,f2,b2,c2,f2);

			make_inst(FMACPIPE,FMADD,e3,b3,d3,e3);
			make_inst(FMACPIPE,FNMSUB,f3,b3,c3,f3);
		}
	} else {
		queue_conjmadd(e1,f1,a1,b1,c1,d1);
		queue_conjmadd(e2,f2,a2,b2,c2,d2);
		queue_conjmadd(e3,f3,a3,b3,c3,d3);
	}
}


void complex_three_conjmaddtos(
	int h1,int e1,int a1,int c1,
	int h2,int e2,int a2,int c2,
	int h3,int e3,int a3,int c3)
{
	if (have_hummer())
	{
		make_inst(SIMDHUMMERPIPE,VMADD_RR_RI,h1,a1,c1,e1);
		make_inst(SIMDHUMMERPIPE,VMADD_RR_RI,h2,a2,c2,e2);
		make_inst(SIMDHUMMERPIPE,VMADD_RR_RI,h3,a3,c3,e3);

		make_inst(SIMDHUMMERPIPE,VMADD_II_MIR,h1,a1,c1,h1);
		make_inst(SIMDHUMMERPIPE,VMADD_II_MIR,h2,a2,c2,h2);
		make_inst(SIMDHUMMERPIPE,VMADD_II_MIR,h3,a3,c3,h3);
	}
	else if (have_mic())
	{
		// TODO
		complex_conjmaddto(h1, e1, a1, c1);
		complex_conjmaddto(h2, e2, a2, c2);
		complex_conjmaddto(h3, e3, a3, c3);
	}
	else
	{
		complex_conjmaddto(h1,e1,a1,c1);
		complex_conjmaddto(h2,e2,a2,c2);
		complex_conjmaddto(h3,e3,a3,c3);
	}
}

void queue_three_conjmaddtos(int h1,int g1,int e1,int f1,int a1,int b1,int c1,int d1,
	int h2,int g2,int e2,int f2,int a2,int b2,int c2,int d2,
	int h3,int g3,int e3,int f3,int a3,int b3,int c3,int d3)
{
	if (have_fmadd)
	{
		if (have_mic()) {
			make_inst(FMOVPIPE, FMOV, h1, e1);
			make_inst(FMOVPIPE, FMOV, g1, f1);
			make_inst(FMOVPIPE, FMOV, h2, e2);
			make_inst(FMOVPIPE, FMOV, g2, f2);
			make_inst(FMOVPIPE, FMOV, h3, e3);
			make_inst(FMOVPIPE, FMOV, g3, f3);

			make_inst(FMACPIPE, FMADD231, h1, a1, c1);
			make_inst(FMACPIPE, FMADD231, g1, a1, d1);
			make_inst(FMACPIPE, FMADD231, h2, a2, c2);
			make_inst(FMACPIPE, FMADD231, g2, a2, d2);
			make_inst(FMACPIPE, FMADD231, h3, a3, c3);
			make_inst(FMACPIPE, FMADD231, g3, a3, d3);

			make_inst(FMACPIPE, FMADD231,  h1, b1, d1);
			make_inst(FMACPIPE, FNMADD231, g1, b1, c1);
			make_inst(FMACPIPE, FMADD231,  h2, b2, d2);
			make_inst(FMACPIPE, FNMADD231, g2, b2, c2);
			make_inst(FMACPIPE, FMADD231,  h3, b3, d3);
			make_inst(FMACPIPE, FNMADD231, g3, b3, c3);
		}
		else {
			make_inst(FMACPIPE,FMADD,h1,a1,c1,e1);
			make_inst(FMACPIPE,FMADD,g1,a1,d1,f1);

			make_inst(FMACPIPE,FMADD,h2,a2,c2,e2);
			make_inst(FMACPIPE,FMADD,g2,a2,d2,f2);

			make_inst(FMACPIPE,FMADD,h3,a3,c3,e3);
			make_inst(FMACPIPE,FMADD,g3,a3,d3,f3);

			make_inst(FMACPIPE,FMADD,h1,b1,d1,h1);
			make_inst(FMACPIPE,FNMSUB,g1,b1,c1,g1);

			make_inst(FMACPIPE,FMADD,h2,b2,d2,h2);
			make_inst(FMACPIPE,FNMSUB,g2,b2,c2,g2);

			make_inst(FMACPIPE,FMADD,h3,b3,d3,h3);
			make_inst(FMACPIPE,FNMSUB,g3,b3,c3,g3);
		}
	}
	else
	{
		queue_conjmaddto(h1,g1,e1,f1,a1,b1,c1,d1);
		queue_conjmaddto(h2,g2,e2,f2,a2,b2,c2,d2);
		queue_conjmaddto(h3,g3,e3,f3,a3,b3,c3,d3);
	}
}

void complex_three_conjmuls(int e1,int a1,int c1,
	int e2,int a2,int c2,
	int e3,int a3,int c3)
{
	if (have_hummer())
	{
		make_inst(SIMDHUMMERPIPE,VMUL_RR_RI,e1,a1,c1);
		make_inst(SIMDHUMMERPIPE,VMUL_RR_RI,e2,a2,c2);
		make_inst(SIMDHUMMERPIPE,VMUL_RR_RI,e3,a3,c3);

		make_inst(SIMDHUMMERPIPE,VMADD_II_MIR,e1,a1,c1,e1);
		make_inst(SIMDHUMMERPIPE,VMADD_II_MIR,e2,a2,c2,e2);
		make_inst(SIMDHUMMERPIPE,VMADD_II_MIR,e3,a3,c3,e3);
	}
	else if (have_mic())
	{
		// TODO
		complex_conjmul(e1, a1, c1);
		complex_conjmul(e2, a2, c2);
		complex_conjmul(e3, a3, c3);
	}
	else
	{
		complex_conjmul(e1,a1,c1);
		complex_conjmul(e2,a2,c2);
		complex_conjmul(e3,a3,c3);
	}
}

void queue_three_conjmuls(int e1,int f1,int a1,int b1,int c1,int d1,
	int e2,int f2,int a2,int b2,int c2,int d2,
	int e3,int f3,int a3,int b3,int c3,int d3)
{

	if (have_fmadd) {
		if (have_mic()) {
			make_inst(FMULPIPE,FMUL,e1,a1,c1);
			make_inst(FMULPIPE,FMUL,f1,a1,d1);
			make_inst(FMULPIPE,FMUL,e2,a2,c2);
			make_inst(FMULPIPE,FMUL,f2,a2,d2);
			make_inst(FMULPIPE,FMUL,e3,a3,c3);
			make_inst(FMULPIPE,FMUL,f3,a3,d3);

			make_inst(FMACPIPE, FMADD231,  e1, b1, d1);
			make_inst(FMACPIPE, FNMADD231, f1, b1, c1);
			make_inst(FMACPIPE, FMADD231,  e2, b2, d2);
			make_inst(FMACPIPE, FNMADD231, f2, b2, c2);
			make_inst(FMACPIPE, FMADD231,  e3, b3, d3);
			make_inst(FMACPIPE, FNMADD231, f3, b3, c3);
		}
		else {
			make_inst(FMULPIPE,FMUL,e1,a1,c1);
			make_inst(FMULPIPE,FMUL,f1,a1,d1);

			make_inst(FMULPIPE,FMUL,e2,a2,c2);
			make_inst(FMULPIPE,FMUL,f2,a2,d2);

			make_inst(FMULPIPE,FMUL,e3,a3,c3);
			make_inst(FMULPIPE,FMUL,f3,a3,d3);

			make_inst(FMACPIPE,FMADD,e1,b1,d1,e1);
			make_inst(FMACPIPE,FNMSUB,f1,b1,c1,f1);

			make_inst(FMACPIPE,FMADD,e2,b2,d2,e2);
			make_inst(FMACPIPE,FNMSUB,f2,b2,c2,f2);

			make_inst(FMACPIPE,FMADD,e3,b3,d3,e3);
			make_inst(FMACPIPE,FNMSUB,f3,b3,c3,f3);
		}
	}
	else
	{
		queue_conjmul(e1,f1,a1,b1,c1,d1);
		queue_conjmul(e2,f2,a2,b2,c2,d2);
		queue_conjmul(e3,f3,a3,b3,c3,d3);
	}
}


/**
 * d = (a*b)+d
 */
void queue_fmac(int d, int a, int b)
{
	queue_fmadd(d, a, b, d);
}


void queue_fnmac(int d, int a, int b)
{
	queue_fnmsub(d, a, b, d);
}


/**
 * scalar floating-point MSUB: d = (a*b) - c.
 */
void queue_fmsub(int d, int a, int b, int c)
{
	int t;

	if (have_fmadd) {
		if (have_mic()) {
			// must reduce 4=>3 operands for MIC
			if (a == d)      make_inst(FMACPIPE, FMSUB213, d, b, c);
			else if (b == d) make_inst(FMACPIPE, FMSUB213, d, a, c);
			else if (c == d) make_inst(FMACPIPE, FMSUB231, d, a, b);
			else {
				make_inst(FMOVPIPE, FMOV, d, a);
				make_inst(FMACPIPE, FMSUB213, d, b, c);
			}
		}
		else {
			make_inst(FMACPIPE, FMSUB, d, a, b, c);
		}
	}
	else {
		t = get_rotating_register(CMADregs);
		make_inst(FMULPIPE, FMUL, t, a, b);
		make_inst(FADDPIPE, FSUB, d, t, c);
	}
}


/**
 * scalar floating-point NMSUB: d = c - (a*b).
 */
void queue_fnmsub(int d, int a, int b, int c)
{
	int t;

	if (have_fmadd)
	{
		if (have_mic()) {
			// must reduce 4=>3 operands for MIC
			if (a == d)      make_inst(FMACPIPE, FNMADD213, d, b, c);
			else if (b == d) make_inst(FMACPIPE, FNMADD213, d, a, c);
			else if (c == d) make_inst(FMACPIPE, FNMADD231, d, a, b);
			else {
				make_inst(FMOVPIPE, FMOV, d, a);
				make_inst(FMACPIPE, FNMADD213, d, b, c);
			}
		}
		else {
			make_inst(FMACPIPE, FNMSUB, d, a, b, c);
		}
	}
	else
	{
		t = get_rotating_register(CMADregs);
		make_inst(FMULPIPE, FMUL, t, a, b);
		make_inst(FADDPIPE, FSUB, d, c, t);
	}
}


/**
 * Scalar floating-point MADD: d = (a*b)+c
 */
void queue_fmadd(int d, int a, int b, int c)
{
	int t;

	if (have_fmadd)
	{
		if (have_mic()) {
			// must reduce 4 operands to 3 operands
			if (d == a)      make_inst(FMACPIPE, FMADD213, d, b, c);
			else if (d == b) make_inst(FMACPIPE, FMADD213, d, a, c);
			else if (d == c) make_inst(FMACPIPE, FMADD231, d, a, b);
			else {
				make_inst(FMOVPIPE, FMOV, d, a);
				make_inst(FMACPIPE, FMADD213, d, b, c);
			}
		}
		else make_inst(FMACPIPE, FMADD, d, a, b, c);
	}
	else
	{
		t = get_rotating_register(CMADregs);
		make_inst(FMULPIPE, FMUL, t, a, b);
		make_inst(FADDPIPE, FADD, d, c, t);
	}
}


/*Bugger is there are no immediate offsets on the damn hummer*/
/*this forces endless loading*/
/**
 * Loads a vector complex from Himm(basereg) into \p dest.
 */
void complex_load(int dest, int Himm, int basereg, enum Datum t)
{
	if (have_hummer())
	{
	  int index_reg;
	  index_reg = pregs[get_constant(get_offset(Himm))];
	  if (t == Double)
	    make_inst(FLODPIPE,CLOAD_INDEXED,dest,index_reg,basereg);
	  else
	    make_inst(FLODPIPE,CSLOAD_INDEXED,dest,index_reg,basereg);
	}
	else if (have_mic())
	{
		if (t == Double)
			make_inst(FLODPIPE, CLOAD_INDEXED, dest, Himm, basereg);
		else
			make_inst(FLODPIPE, CSLOAD_INDEXED, dest, Himm, basereg);
	}
	else
	{
		for (int v = 0; v < nsimd()*2; v++)
		{
			int destv = dest+v;
			int units = get_offset(Himm)+ v*SizeofDatum(t);
			int Himmv = get_offset_handle(units, Byte);
			queue_fload(destv,Himmv,basereg,t);
		}
	}
}

void complex_rsplat(int dest, int Himm, int basereg, enum Datum t)
{
  int index_reg;
  if (have_hummer()){

    queue_fload(dest,Himm,basereg,t);
    int imm = def_offset(0,Byte,"BGQ_RSPLAT");
    make_inst(FMOVPIPE,FSPLAT,dest,dest,imm);

  } else if (have_mic()) {

    queue_fload(dest,Himm,basereg,t); // already does a splat on knc

  } else {

    for (int v = 0; v < nsimd()*2; v++) {
      int destv = dest+v;
      queue_fload(destv,Himm,basereg,t);
    }

  }
}

void complex_csplat(int dest, int Himm, int basereg, enum Datum t)
{
  int index_reg;
  if (have_hummer()){

    int index_reg;
    index_reg = pregs[get_constant(get_offset(Himm))];
    make_inst(FLODPIPE,CLOAD_SPLAT_INDEXED,dest,index_reg,basereg); // does splat on bgq

  } else if (have_mic()) {

    // probably Bad code -- what if Himm varies? Code will fail with collision in names.
    if ( t==Double){
      int o=get_offset(Himm);
      o=o+8;
      int Himmpp = def_offset(8,Byte,"cpslat_8");
      make_inst(SIMDMICSWIZPIPE,CLOAD_SPLAT,dest,MIC_MASK_ODD,Himm,basereg);  // does splat on real
      make_inst(SIMDMICSWIZPIPE,CLOAD_SPLAT,dest,MIC_MASK_EVEN,Himmpp,basereg); // does splat on imag
    } else {
      int o=get_offset(Himm);
      o=o+4;
      int Himmpps = def_offset(8,Byte,"cpslat_4");
      make_inst(SIMDMICSWIZPIPE,CSLOAD_SPLAT,dest,MIC_MASK_ODD,Himm,basereg);  // does splat on real
      make_inst(SIMDMICSWIZPIPE,CSLOAD_SPLAT,dest,MIC_MASK_EVEN,Himmpps,basereg); // does splat on imag
    }

  } else {

    for (int v = 0; v < nsimd()*2; v++) {
      int destv = dest+v;
      queue_fload(destv,Himm,basereg,t);
    }

  }
}



/**
 * Stores vector complex: src => Himm(basereg)
 */
void complex_store(int src, int Himm, int basereg, enum Datum t)
{
	int index_reg;

	if (have_hummer())
	{
		index_reg = pregs[get_constant(get_offset(Himm))];
		if (t == Double) {
		  make_inst(FSTRPIPE,CSAVE_INDEXED,src,index_reg,basereg);
		} else if ( t == Single ) {
		  make_inst(FSTRPIPE,CSSAVE_INDEXED,src,index_reg,basereg);
		}else {
		  exit(-1);
		}
	}
	else if (have_mic())
	{
		if (t == Double)
			make_inst(FSTRPIPE, CSAVE_INDEXED, src, Himm, basereg);
		else
			make_inst(FSTRPIPE, CSSAVE_INDEXED, src, Himm, basereg);
	}
	else
	{
		for (int v = 0; v < nsimd()*2; v++)
		{
			int srcv = src + v;
			int units = get_offset(Himm) + v*SizeofDatum(t);
			int Himmv = get_offset_handle(units, Byte);
			queue_fstore(srcv, Himmv, basereg, t);
		}
	}
}

/*
 */
void complex_load_half(int dest,int Himm,int basereg,int membuf,int t1,int t2, int mask)
{
  int zero_reg = pregs[get_constant(0)];
  queue_iload(t1,Himm,basereg);
  queue_andc(t2,t1,mask);
  queue_and(t1,t1,mask);
  queue_lshift(t2,t2,get_offset_handle(16,Byte));
  queue_istore(t1,get_offset_handle(0,Byte),membuf);       // Load 128 bits
  queue_istore(t2,get_offset_handle(8,Byte),membuf);
  make_inst(DIRECTIVE,LS_BARRIER);
  make_inst(FLODPIPE,CSLOAD_INDEXED,dest,zero_reg,membuf);// Save 128 bits
}
void complex_store_half(int src ,int Himm,int basereg,int membuf,int t1,int t2, int mask)
{
  int zero_reg = pregs[get_constant(0)];
  make_inst(FSTRPIPE,CSSAVE_INDEXED,src,zero_reg,membuf);// Save 128 bits
  make_inst(DIRECTIVE,LS_BARRIER);
  queue_iload(t1,get_offset_handle(0,Byte),membuf);       // Load 128 bits
  queue_iload(t2,get_offset_handle(8,Byte),membuf);
  queue_and(t1,t1,mask);
  queue_and(t2,t2,mask);
  queue_rshift(t2,t2,get_offset_handle(16,Byte));
  make_inst(IALUPIPE, IOR, t1, t2, t1);
  queue_istore(t1,Himm,basereg);
}
// QUESTION: what's the diff between complex_add() and simd_add()?
/**
 * Vector complex ADD: e = a + c
 */
void complex_add(int e, int a, int c)
{
	if (have_hummer())
	{
		make_inst(SIMDHUMMERPIPE, VADD, e, a, c);
	}
	else if (have_mic())
	{
		make_inst(SIMDMICPIPE, MIC_ADD, e, a, c);
	}
	else
	{
		for (int v = 0; v < 2*nsimd(); v++)
		{
			queue_fadd(e+v, a+v, c+v);
		}
	}
}


/**
 * Vector complex SUBTRACT: e = a - c
 */
void complex_sub(int e, int a, int c)
{
	if (have_hummer())
	{
		make_inst(SIMDHUMMERPIPE, VSUB, e, a, c);
	}
	else if (have_mic())
	{
		make_inst(SIMDMICPIPE, MIC_SUB, e, a, c);
	}
	else
	{
		for (int v = 0; v < 2*nsimd(); v++)
		{
			queue_fsub(e+v, a+v, c+v);
		}
	}
}


/**
 * (e+if) = (a+ib) + i*(c+id)
 */
void complex_ApiB(int e, int a, int c)
{
	if (have_hummer())
	{
		complex_constants_assert(cmplx_i);
		make_inst(SIMDHUMMERPIPE, VMADD_MII_IR, e, reg_tmp, c, a);
	}
	else if (have_mic())
	{
		make_inst(SIMDMICSWIZPIPE, MIC_SUB_MASK_SWIZ, e, MIC_MASK_ODD,  a, c, MIC_SWIZ_CDAB);  // e = a-d
		make_inst(SIMDMICSWIZPIPE, MIC_ADD_MASK_SWIZ, e, MIC_MASK_EVEN, a, c, MIC_SWIZ_CDAB);  // f = b+c
	}
	else
	{
		for(int v = 0; v < nsimd(); v++)
		{
			int b = a+1, f = e+1, d = c+1;
			queue_fsub(e, a, d);
			queue_fadd(f, b, c);
			a+=2; e+=2; c+=2;
		}
	}
}


/**
 * (e+if) = (a+ib) - i(c+id)
 *        = (a+d)  + i(b-c)
 */
void complex_AmiB(int e, int a, int c)
{
	if (have_hummer())
	{
		if (PROC->cross_negative)
		{
			complex_constants_assert(cmplx_i);
			make_inst(SIMDHUMMERPIPE,VMADD_II_MIR,e,reg_tmp,c,a);
//      make_inst(SIMDHUMMERPIPE,VMADD_II_MIR,e,c,reg_tmp,a);
		}
		else
		{
			complex_constants_assert(cmplx_neg_i);
			make_inst(SIMDHUMMERPIPE,VMADD_MII_IR,e,reg_tmp,c,a);
		}
	}
	else if (have_mic())
	{
		make_inst(SIMDMICSWIZPIPE, MIC_ADD_MASK_SWIZ, e, MIC_MASK_ODD,  a, c, MIC_SWIZ_CDAB);  // e = a+d
		make_inst(SIMDMICSWIZPIPE, MIC_SUB_MASK_SWIZ, e, MIC_MASK_EVEN, a, c, MIC_SWIZ_CDAB);  // f = b-c

		/*
		alternative: this needs a copy and a sign vector reg.
		result = (a,b) + (S_cdab(c,d) * (+1,-1))
		SIGN_PM is a vector reg with (+1,-1,+1,-1...)
		One input operand must be clobbered; this will be the sign which is copied beforehand.

		make_inst(SIMDMICPIPE, MIC_MOV, e, SIGN_PM);  // copy sign (possibly from mem?)
		make_inst(SIMDMICSWIZPIPE, MIC_FMADD132_SWIZ, e, a, c, MIC_SWIZ_CDAB);
		*/
	}
	else
	{
		for (int v = 0; v < nsimd(); v++)
		{
			int b = a+1, f = e+1, d = c+1;
			queue_fadd(e, a, d);
			queue_fsub(f, b, c);
			a+=2; e+=2; c+=2;
		}
	}
}
/* ...END complex queueing routines*/


/*
* Provided as general primitives for treating SIMD index
* as a 2^n grid. I need to clarify my thoughts to code this
* for general vector length.
*/
// permute patterns
static int const PATT_PERM = (0x2<<9)|(0x3<<6)|(0x0<<3)|0x1;  // 2,3,0,1
static int const PATT_LOLO = (0x0<<9)|(0x1<<6)|(0x4<<3)|0x5;  // 0,1,4,5
static int const PATT_HIHI = (0x2<<9)|(0x3<<6)|(0x6<<3)|0x7;  // 2,3,6,7


// permute handles
static int hperm, hlolo, hhihi; // BGQ
static int mic_perm_2301,mic_perm_1032; //KNC

static int permreg = -1;
static int permval = -1;


/**
 * Creates permute handles and sets permute register.
 */
void complex_simd_init(int preg)
{
	hperm  = def_offset(PATT_PERM, Byte, "perm");
	hlolo  = def_offset(PATT_LOLO, Byte, "lolo");
	hhihi  = def_offset(PATT_HIHI, Byte, "hihi");

	// Intel : 3210 is identity map
	mic_perm_2301 = def_offset(0xB1,Byte,"mic_perm_2301"); //b1011,0001 = B,3,0,1
	mic_perm_1032 = def_offset(0x4E,Byte,"mic_perm_1032"); //b0100,1110 = 4E
	permreg = preg;
}

void complex_simd_permute(int mu, int dest, int src)
{
	if (nsimd()%2 != 0) {
		fprintf(stderr, "Error: SIMD vectors cannot hold 4n elements.\n");
		exit(0);
	}

	if ( have_hummer() ) {
	  if (mu != 3) mu=3;
	  complex_simd_permute_internal(dest, src, src, hperm);
	}
	if ( have_mic() ) {
	  if (mu==3) {
	    fprintf(stderr,"calling vpermf32x4 mu=3\n");
	    //	    make_inst(SIMDMICSWIZPIPE, MIC_VPERMF32X4, dest, MIC_MASK_FULL, src, mic_perm_1032);
	    make_inst(SIMDPERMUTEPIPE, VPERMCOARSE, dest, mic_perm_1032, src);

	  } else if (mu==2) {
	    fprintf(stderr,"calling vpermf32x4 mu=2\n");
	    //	    make_inst(SIMDMICSWIZPIPE, MIC_VPERMF32X4, dest, MIC_MASK_FULL, src, mic_perm_2301);
	    make_inst(SIMDPERMUTEPIPE, VPERMCOARSE, dest, mic_perm_2301, src);
	  } else {
	    fprintf(stderr,"calling blend mu=1\n");
	    make_inst(SIMDMICSWIZPIPE, MIC_BLEND_MASK_SWIZ, dest, MIC_MASK_FULL, src, src, MIC_SWIZ_BADC);
	  }          
	}
}
void complex_simd_blend(int mu,int dest,int src1,int src2)
{
  int mask,swiz;
  if(mu==3){
    mask = MIC_MASK_MU3; // 0xF0
    swiz = MIC_SWIZ_DCBA;
    make_inst(SIMDMICSWIZPIPE, MIC_BLEND_MASK_SWIZ, dest, mask, src1, src2, swiz);
  } else {
    mask = MIC_MASK_MU2; // 0xCC
    swiz = MIC_SWIZ_DCBA; 
    make_inst(SIMDMICSWIZPIPE, MIC_BLEND_MASK_SWIZ, dest, mask, src1, src2, swiz);
  }
}

void complex_simd_extract(int hi, int mu, int dest , int src1 , int src2)
{
	if (nsimd()%2 != 0) {
		fprintf(stderr, "Error: SIMD vectors cannot hold 4n elements.\n");
		exit(0);
	}

	if ( have_hummer() ) {
	  if (mu != 3) exit(0);

	  if (hi) complex_simd_permute_internal(dest, src1, src2, hhihi);
	  else    complex_simd_permute_internal(dest, src1, src2, hlolo);
	}

	if ( have_mic() ) {
	  if (hi) {
	    complex_simd_permute(mu,dest,src1);
	    complex_simd_blend(mu,dest,dest,src2);
	  } else {
	    complex_simd_permute(mu,dest,src2);
	    complex_simd_blend(mu,dest,dest,src1);
	  }
	}
}


void complex_simd_merge(int hi, int mu, int dest1, int src1, int src2)
{
  complex_simd_extract(hi,mu,dest1,src1,src2);
}


/**
 * \param hpatt Handle of permute pattern.
 */
void complex_simd_permute_internal(int dest, int src1, int src2, int hpatt)
{
	int i, pi, tmp;
	int perm = get_offset(hpatt);
	if (have_permute()) {

		if (have_mic()) {
			// bloody MIC hack. We cannot directly use the hpatt
			int swiz, mask;
			switch (perm) {
				case PATT_PERM:
					mask = MIC_MASK_FULL;
					swiz = MIC_SWIZ_BADC;
					break;
				case PATT_LOLO:
					mask = MIC_MASK_EVEN_PAIRS;
					swiz = MIC_SWIZ_BADC;
					break;
				case PATT_HIHI:
					tmp = src1;  // NB: must swap sources!
					src1 = src2;
					src2 = tmp;
					mask = MIC_MASK_ODD_PAIRS;
					swiz = MIC_SWIZ_BADC;
					break;
				default:
					fprintf(stderr, "Error: Invalid permute pattern.\n");
					abort();
			}
			make_inst(SIMDMICSWIZPIPE, MIC_BLEND_MASK_SWIZ, dest, mask, src1, src2, swiz);
		}
		else {
			if (permreg == -1) {
				fprintf(stderr, "Error: Permute reg not set.\n");
				abort();
			}

			if (permval != perm) {
				make_inst(SIMDPERMUTEPIPE,VPERMIMM,permreg,hpatt);
				permval = perm;
			}
			make_inst(SIMDPERMUTEPIPE,VPERM,dest,src1,src2,permreg);
		}
	}
	else
	{
		if (permreg == -1) {
			fprintf(stderr, "Error: Permute reg not set.\n");
			abort();
		}

		for (i=0; i<2*nsimd(); i++) {
			// group vector elements in units of 4
			pi = (perm>>(3*(3-(i%4))))&0x7;  // get 3 LSBs after right-shift
			if (pi < 4) queue_fmov(permreg+i,src1+pi+(i/4));
			else        queue_fmov(permreg+i,src2+pi-4+(i/4));
		}
		for (i=0; i<2*nsimd(); i++) {
			queue_fmov(dest+i,permreg+i);
		}
	}
}
void complex_zero(int reg)
{
  if (have_hummer()) {
    make_inst(SIMDHUMMERPIPE,VZERO,reg);
  } else { 
    complex_load(reg,imm_0,i_addr,Double);
  }
}
