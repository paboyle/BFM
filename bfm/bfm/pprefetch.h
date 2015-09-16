/*------------------------------------------------------------------*/
/*                                                                  */
/* (C) Copyright IBM Corp.  2009, 2009                              */
/* Eclipse Public License (EPL)                                     */
/*                                                                  */
/*------------------------------------------------------------------*/
#ifndef _L1P_PERFECT_H_ // Prevent multiple inclusion.
#define _L1P_PERFECT_H_

#include <hwi/include/common/compiler_support.h>

__BEGIN_DECLS

#include <sys/types.h>
#include <stdint.h>
#include <stdlib.h>
#include <spi/include/kernel/memory.h>
#include <spi/include/kernel/location.h>
#include <spi/include/l1p/types.h>
#include <hwi/include/bqc/l1p_mmio.h>
#include <hwi/include/bqc/A2_inlines.h>

extern L1P_SPIContext _L1P_Context[68];

#undef DD1_L1P_Workaround_local

/*!
 \brief Allocates enough storage so that the perfect prefetcher can track up to <n> L1 misses.  
 \param[in] n The maximum number of L1 misses that can be tracked by the list.
 
 Storage is retained until:
 1)	L1P_Unconfigure() is performed.
 2)	L1P_SetPattern() is performed.
 
 If  the L1P_Configure() command is nested: 
 *	If nesting mode has been set to L1P_NestingSaveContext then the L1P SPI will push a L1P context structure onto a stack of L1P context structures.  
                When an L1P_Unconfigure() is called, this L1P context structure will be restored.  This is the default mode.
 *	If nesting mode has been set to L1P_NestingIgnore, then the L1P SPI will reference count the L1P_Configures.  When nested, the SPI will not write 
                new pattern addresses into the L1p hardware.  Once the same number of L1P_Unconfigure() routines have been called, the L1P SPI will return to normal function.  
 *	If the nesting node has been set to L1P_NestingFlat, then the L1P SPI will be reference count and ignore the L1P_Configures.  Any L1P_SetPattern() 
                calls will be ignored if they occur in a nested context.
 *	If nesting mode has been set to L1P_NestingError, the L1P SPI will display an error message and assert.  This will terminate the active process 
                with a corefile.  This mode is to be used for debug purposes.  
 */
#define L1P_LISTWORDSIZE (4)
__INLINE__ int L1P_PatternConfigure(uint64_t n)
{
    L1P_SPIContext* context = &_L1P_Context[Kernel_ProcessorID()];
    
    // \todo Add support for nesting L1P_PatternConfigures.
    
    context->currentPattern.size = n*L1P_LISTWORDSIZE;
    context->currentPattern.ReadPattern = malloc(context->currentPattern.size+128);
    if(context->currentPattern.ReadPattern == NULL)
        return -1;
    memset(context->currentPattern.ReadPattern,255,context->currentPattern.size);
    context->currentPattern.WritePattern = malloc(context->currentPattern.size+128);
    if(context->currentPattern.WritePattern == NULL)
        return -1;   
    memset(context->currentPattern.WritePattern,255,context->currentPattern.size);
    context->implicitPatternAllocate = 1;
    
    *((uint64_t*)(Kernel_L1pBaseAddress() + L1P_PP_THREADOFFSET * Kernel_ProcessorThreadID() + L1P_PP_CTRL))   = L1P_PP_CTRL_DEPTH(7) | L1P_PP_CTRL_MAXLIST(7);
    *((uint64_t*)(Kernel_L1pBaseAddress() + L1P_PP_THREADOFFSET * Kernel_ProcessorThreadID() + L1P_PP_MAXTOL)) = 0xffffffff;
    
    return Kernel_L1pSetPatternAddress(context->currentPattern.ReadPattern, context->currentPattern.WritePattern, context->currentPattern.size);
}

/*!
 \brief Deallocates storage used by the L1p SPI.
 If one is available, the L1P SPI will pop a L1P context structure from the stack of L1P context structures.  The context will then be used to restore the previous L1P pattern status and pointers.  
 */
__INLINE__ int L1P_PatternUnconfigure()
{
    int rc;
    L1P_SPIContext* context = &_L1P_Context[Kernel_ProcessorID()];
    // \todo Add support for nesting L1P_PatternConfigures.
    // \todo Add error checking.  Unconfigure w/o configure?

    if(context->implicitPatternAllocate)
    {
        rc = Kernel_L1pSetPatternAddress(NULL, NULL, 0);
        if(rc) return rc;
        
        if(context->currentPattern.ReadPattern == NULL)
            free(context->currentPattern.ReadPattern);
        if(context->currentPattern.WritePattern == NULL)
            free(context->currentPattern.WritePattern);
        context->currentPattern.ReadPattern = NULL;
        context->currentPattern.WritePattern = NULL;
        context->currentPattern.size = 0;
    }
    return 0;
}

/*!
 \brief The perfect prefetcher will start monitoring L1 misses and performing prefetch requests based on those misses.  
 
 \param[in] record Instructs the PatternStart to record the pattern of L1 misses for the next iteration.  
 This L1P_PatternStart() should be called at the beginning of every entrance into the section of code that has been recorded.  
 */
__INLINE__ int L1P_PatternStart(int record)
{
  *((uint64_t*)(Kernel_L1pBaseAddress() + L1P_PP_THREADOFFSET * Kernel_ProcessorThreadID() + L1P_PP_STATUS)) = 0xf;
  *((uint64_t*)(Kernel_L1pBaseAddress() + L1P_PP_THREADOFFSET * Kernel_ProcessorThreadID() + L1P_PP_CTRL)) |= L1P_PP_CTRL_START | L1P_PP_CTRL_LOAD | ((record == 0)?L1P_PP_CTRL_INHIBIT:0);
  return 0;
}

/*!
 \brief Suspends the active perfect prefetcher.  
 The Linear Stream Prefetcher and the other 3 Perfect Prefetchers on the core will continue to execute.
 This routine can be used in conjunction with L1P_PatternResume() to avoid recording out-of-bound memory fetches - - such as instructions performing a periodic printf.    
 It can also be used to avoid sections of code that perform memory accesses that are inconsistent between iterations.  
 */
__INLINE__ int L1P_PatternPause()
{
    *((uint64_t*)(Kernel_L1pBaseAddress() + L1P_PP_THREADOFFSET * Kernel_ProcessorThreadID() + L1P_PP_CTRL)) |= L1P_PP_CTRL_PAUSE;
    return 0;
}

/*! 
 \brief Resumes the perfect prefetcher from the last pattern offset location.  
 This routine can be used in conjunction with L1P_PatternPause() to avoid recording memory fetches that are not likely to repeat - - such as instructions 
 performing a periodic printf.    It can also be used to avoid sections of code that perform memory accesses that are inconsistent between iterations
*/
__INLINE__ int L1P_PatternResume()
{
    *((uint64_t*)(Kernel_L1pBaseAddress() + L1P_PP_THREADOFFSET * Kernel_ProcessorThreadID() + L1P_PP_CTRL)) &= ~L1P_PP_CTRL_PAUSE;
    return 0;
}

/*!
 \brief Stops the perfect prefetcher and resets the list offsets to zero
 */
__INLINE__ int L1P_PatternStop()
{    
    uint64_t ctrl = *((uint64_t*)(Kernel_L1pBaseAddress() + L1P_PP_THREADOFFSET * Kernel_ProcessorThreadID() + L1P_PP_CTRL));
#if DD1_L1P_Workaround_local
    int x;
    int offset = Kernel_ProcessorCoreID()*4;
    uint32_t freevalue = 0;
    uint64_t my_unique_index = Kernel_ProcessorID() + 1;
    for(x=0; x<4; x++)
    {
        while(1)
        {
            if(Compare_and_Swap32(&_L1P_Context[x + offset].criticalatom, &freevalue, my_unique_index))
            {
                break;
            }
            ThreadPriority_Low();
            while (_L1P_Context[x + offset].criticalatom) { }
            ThreadPriority_Medium(); // Use high priority while we are holding the kernel lock                                                                
            freevalue = 0;
        }
        *((uint64_t*)(Kernel_L1pBaseAddress() + L1P_PP_THREADOFFSET * x + L1P_PP_CTRL)) |= L1P_PP_CTRL_PAUSE;
    }
    ppc_msync();
#else 
    static int suppress;
    if ( !suppress ) printf("not using DD1 workaround\n");
    suppress=1;
#endif
    *((uint64_t*)(Kernel_L1pBaseAddress() + L1P_PP_THREADOFFSET * Kernel_ProcessorThreadID() + L1P_PP_CTRL)) = ctrl & ~(L1P_PP_CTRL_START | L1P_PP_CTRL_LOAD | L1P_PP_CTRL_INHIBIT);
    
    // poll the list
    // finish & abandon & max_length_reach & end_of_list
    uint64_t list_status=0;
    while( (list_status & L1P_PP_STATUS_FINISHED) == 0 ) 
    {
        list_status = in64((uint64_t*)(Kernel_L1pBaseAddress() + L1P_PP_THREADOFFSET * Kernel_ProcessorThreadID() +  L1P_PP_STATUS));
    }
#if DD1_L1P_Workaround_local
    for(x=0; x<4; x++)
    {
        *((uint64_t*)(Kernel_L1pBaseAddress() + L1P_PP_THREADOFFSET * x + L1P_PP_CTRL)) &= (~L1P_PP_CTRL_PAUSE);
        ppc_msync();
        
        // Consider checking that we own this lock before resetting it. This test may not be necessary if we can trust ourself!                                 
        // reset the lock                                                                                                                                       
        _L1P_Context[x + offset].criticalatom = 0;
        // Move our priority back down to normal. We were running with a high priority while holding the kernel lock                                            
        ThreadPriority_Medium(); // Use high priority while we are holding the kernel lock                                                                      
    }
    ppc_msync();
#endif

    if((ctrl & L1P_PP_CTRL_INHIBIT) == 0)
    {
        L1P_SPIContext* context = &_L1P_Context[Kernel_ProcessorID()];
        // swap generate and record lists to setup for next iteration
        
        void* tmp = context->currentPattern.ReadPattern;
        context->currentPattern.ReadPattern = context->currentPattern.WritePattern;
        context->currentPattern.WritePattern = tmp;
        Kernel_L1pSetPatternAddress(context->currentPattern.ReadPattern, context->currentPattern.WritePattern, context->currentPattern.size);
    }
    return 0;
}

/*!
 \brief Returns the current status for the L1 perfect prefetcher
 \param[out] status The current L1p status will be placed at this memory location.
 */
__INLINE__ int L1P_PatternStatus(L1P_Status_t* status)
{
    uint64_t tmp = *((uint64_t*)(Kernel_L1pBaseAddress() + L1P_PP_THREADOFFSET * Kernel_ProcessorThreadID() + L1P_PP_STATUS));
    status->status = tmp;
    return 0;
}

/*! 
 \brief Returns the current pattern depths for the L1 perfect prefetcher.  
 The pattern depth is the current index into the pattern that the L1p is executing.  
 \param[out] fetch_depth Used to determine how far in the current pattern/sequence the L1p has progressed.  
 \param[out] generate_depth Used to optimize the pattern length parameter to L1P_PatternConfigure() in order to reduce the memory footprint of the L1p pattern
 */
__INLINE__ int L1P_PatternGetCurrentDepth(uint64_t* fetch_depth, uint64_t* generate_depth)
{
    int rc;
    L1P_SPIContext* context = &_L1P_Context[Kernel_ProcessorID()];
    void* curread;
    void* curwrite;
    rc = Kernel_L1pGetCurrentPatternAddress(&curread, &curwrite);
    if(rc) return rc;
    *fetch_depth    = (uint64_t)curread - (uint64_t)context->currentPattern.ReadPattern;
    *generate_depth = (uint64_t)curwrite - (uint64_t)context->currentPattern.WritePattern;
    return 0;
}

/*!
 \brief Returns the current nesting mode for the L1 perfect prefetcher.
 \see L1P_PatternNest_t
 \see L1P_PatternSetNestingMode
 */
__INLINE__ int L1P_PatternGetNestingMode(L1P_PatternNest_t* mode)
{
    L1P_SPIContext* context = &_L1P_Context[Kernel_ProcessorID()];
    *mode = context->nestPolicy;
    return 0;
}

/*!
 \brief Sets the nesting mode for the L1 perfect prefetcher.
 \see L1P_PatternNest_t
 \see L1P_PatternSetNestingMode
 */
__INLINE__ int L1P_PatternSetNestingMode(L1P_PatternNest_t mode)
{
    // \todo validate mode
    L1P_SPIContext* context = &_L1P_Context[Kernel_ProcessorID()];
    context->nestPolicy = mode;
    return 0;
}

/*!
 \brief Sets the number of consecutive L1 misses that did not match the current location in the pattern.  
 Once this number has been exceeded, the prefetching activity will cease and the pattern will be marked as "Abandoned" in the L1P_Status_t structure returned by L1P_PatternStatus().
 */
__INLINE__ int L1P_PatternSetAbandonThreshold(uint64_t numL1misses)
{
    *((uint64_t*)(Kernel_L1pBaseAddress() + L1P_PP_THREADOFFSET * Kernel_ProcessorThreadID() + L1P_PP_MAXTOL)) = numL1misses;
    return 0;
}

/*!
 \brief Returns the number of consecutive L1 misses that did not match the current location in the pattern.  
 Once this number has been exceeded, the prefetching activity will cease and the pattern will be marked as "Abandoned" in the L1P_Status_t structure returned by L1P_PatternStatus().  
 */
__INLINE__ int L1P_PatternGetAbandonThreshold(uint64_t* numL1misses)
{
    *numL1misses = *((uint64_t*)(Kernel_L1pBaseAddress() + L1P_PP_THREADOFFSET * Kernel_ProcessorThreadID() + L1P_PP_MAXTOL));
    return 0;
}
 
/*!
 \brief Sets a software enable/disable for L1p perfect prefetecher.  
 This can be used to ascertain whether the usage of the prefetcher is improving performance
 */
__INLINE__ int L1P_PatternSetEnable(int enable)
{
    // \todo Implement L1P_PatternSetEnable
    return -1;
}

/*!
 \brief Returns the software enable/disable for L1p perfect prefetcher
 */
__INLINE__ int L1P_PatternGetEnable(int* enable)
{
    // \todo Implement L1P_PatternGetEnable
    return -1;
}

__END_DECLS

#endif // Add nothing below this line.
