#include <stdio.h>
#include <stdlib.h>
#include <qmp.h>
extern "C" { 
   void* memalign(size_t, size_t);
};

#include <bfm.h>
#include <bfmcommqmp.h>
#include <bfmcommspi.h>

//////////////////////////////////////////////////////////
// SPI thingy-me-bob's
// Save a 32MB region for our comms buffers??
//////////////////////////////////////////////////////////


#include <spi/include/mu/Descriptor.h>
#include <spi/include/mu/Descriptor_inlines.h>
#include <spi/include/mu/InjFifo.h>
#include <spi/include/mu/Addressing.h>
#include <spi/include/mu/Addressing_inlines.h>
#include <spi/include/mu/GIBarrier.h>
#include <spi/include/kernel/MU.h>
#include <hwi/include/bqc/classroute.h>
#include <spi/include/kernel/collective.h>
#include <firmware/include/personality.h>


/////////////////////////////////////////////////////
// Torus mapping related
/////////////////////////////////////////////////////
static uint64_t BfmSPI_Torus_FIFO_Map[10] = {
  MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_AM,
  MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_AP,
  MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_BM,
  MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_BP,
  MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_CM,
  MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_CP,
  MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_DM,
  MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_DP,
  MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_EM,
  MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_EP
};

static uint8_t BfmSPI_hints_abcd[10] = {
  MUHWI_PACKET_HINT_AM|MUHWI_PACKET_HINT_B_NONE|MUHWI_PACKET_HINT_C_NONE|MUHWI_PACKET_HINT_D_NONE,
  MUHWI_PACKET_HINT_AP|MUHWI_PACKET_HINT_B_NONE|MUHWI_PACKET_HINT_C_NONE|MUHWI_PACKET_HINT_D_NONE,
  MUHWI_PACKET_HINT_A_NONE|MUHWI_PACKET_HINT_BM|MUHWI_PACKET_HINT_C_NONE|MUHWI_PACKET_HINT_D_NONE,
  MUHWI_PACKET_HINT_A_NONE|MUHWI_PACKET_HINT_BP|MUHWI_PACKET_HINT_C_NONE|MUHWI_PACKET_HINT_D_NONE,
  MUHWI_PACKET_HINT_A_NONE|MUHWI_PACKET_HINT_B_NONE|MUHWI_PACKET_HINT_CM|MUHWI_PACKET_HINT_D_NONE,
  MUHWI_PACKET_HINT_A_NONE|MUHWI_PACKET_HINT_B_NONE|MUHWI_PACKET_HINT_CP|MUHWI_PACKET_HINT_D_NONE,
  MUHWI_PACKET_HINT_A_NONE|MUHWI_PACKET_HINT_B_NONE|MUHWI_PACKET_HINT_C_NONE|MUHWI_PACKET_HINT_DM,
  MUHWI_PACKET_HINT_A_NONE|MUHWI_PACKET_HINT_B_NONE|MUHWI_PACKET_HINT_C_NONE|MUHWI_PACKET_HINT_DP
};

static uint8_t BfmSPI_hints_e[10] = {
  MUHWI_PACKET_HINT_E_NONE|MUHWI_PACKET_DO_NOT_ROUTE_TO_IO_NODE|MUHWI_PACKET_USE_DETERMINISTIC_ROUTING| MUHWI_PACKET_DO_NOT_DEPOSIT,
  MUHWI_PACKET_HINT_E_NONE|MUHWI_PACKET_DO_NOT_ROUTE_TO_IO_NODE|MUHWI_PACKET_USE_DETERMINISTIC_ROUTING| MUHWI_PACKET_DO_NOT_DEPOSIT,  
  MUHWI_PACKET_HINT_E_NONE|MUHWI_PACKET_DO_NOT_ROUTE_TO_IO_NODE|MUHWI_PACKET_USE_DETERMINISTIC_ROUTING| MUHWI_PACKET_DO_NOT_DEPOSIT,  
  MUHWI_PACKET_HINT_E_NONE|MUHWI_PACKET_DO_NOT_ROUTE_TO_IO_NODE|MUHWI_PACKET_USE_DETERMINISTIC_ROUTING| MUHWI_PACKET_DO_NOT_DEPOSIT,  
  MUHWI_PACKET_HINT_E_NONE|MUHWI_PACKET_DO_NOT_ROUTE_TO_IO_NODE|MUHWI_PACKET_USE_DETERMINISTIC_ROUTING| MUHWI_PACKET_DO_NOT_DEPOSIT,  
  MUHWI_PACKET_HINT_E_NONE|MUHWI_PACKET_DO_NOT_ROUTE_TO_IO_NODE|MUHWI_PACKET_USE_DETERMINISTIC_ROUTING| MUHWI_PACKET_DO_NOT_DEPOSIT,  
  MUHWI_PACKET_HINT_E_NONE|MUHWI_PACKET_DO_NOT_ROUTE_TO_IO_NODE|MUHWI_PACKET_USE_DETERMINISTIC_ROUTING| MUHWI_PACKET_DO_NOT_DEPOSIT,  
  MUHWI_PACKET_HINT_E_NONE|MUHWI_PACKET_DO_NOT_ROUTE_TO_IO_NODE|MUHWI_PACKET_USE_DETERMINISTIC_ROUTING| MUHWI_PACKET_DO_NOT_DEPOSIT,  
  MUHWI_PACKET_HINT_EM    |MUHWI_PACKET_DO_NOT_ROUTE_TO_IO_NODE|MUHWI_PACKET_USE_DETERMINISTIC_ROUTING| MUHWI_PACKET_DO_NOT_DEPOSIT,  
  MUHWI_PACKET_HINT_EP    |MUHWI_PACKET_DO_NOT_ROUTE_TO_IO_NODE|MUHWI_PACKET_USE_DETERMINISTIC_ROUTING| MUHWI_PACKET_DO_NOT_DEPOSIT
};

extern "C" {
extern int                 QMP_Labcde[5];
extern MUHWI_Destination_t QMP_neigh[8];
extern MUHWI_Destination_t QMP_me;
}


static int last_logical_dir;

//////////////////////////////////
// Fast barrier
//////////////////////////////////
static uint32_t BfmSPI_GIBarrierClassRoute = 0; // 0 defaults to comm_world
static MUSPI_GIBarrier_t BfmSPI_GIBarrier;

//////////////////////////////////
// Global reduction
//////////////////////////////////
static MUHWI_Descriptor_t *BfmSPI_gsum_descr;
static int BfmSPI_gsum_buf = 0;


/////////////////////////////////////////////////////////////////////////////
// Memory handling
/////////////////////////////////////////////////////////////////////////////
static uint64_t BFM_DEVICE_MEM_SIZE = (32*1024*1024) ;

#define CHECK(A) { if ( rc ) {printf(A); printf("Looks bad returned %d\n",rc); fflush(stdout); sleep(10) ; exit(rc);}}
//#define CHECK(A) { if ( this->isBoss() ) { printf(A); fflush(stdout);} if ( rc ) { exit(rc);}}

static int BfmSPI_MU_resources_initialised;
static int BfmSPI_MU_resources_reserved;

static Kernel_MemoryRegion_t BfmSPI_device_memory_region;
static unsigned char      *BfmSPI_device_memory;
static unsigned char      *BfmSPI_device_malloc_ptr;
static uint64_t            BfmSPI_device_bytes;

//Restore from these values when we init a BFM object
static unsigned char      *BfmSPI_device_malloc_ptr_save;
static uint64_t            BfmSPI_device_bytes_save;

// Address translation
static uint64_t BfmSPI_BasePa;
static uint64_t BfmSPI_BaseVa;


static  Personality_t personality;
static int isBoss;

// It appears we must spread our allocation thinly over 
// different groups associated with cores. 
static int                              BfmSPI_sub_group[NBRS] = {0,4,8,12,16,20,24,28};
static MUSPI_BaseAddressTableSubGroup_t BfmSPI_BAT_sub_group[NBRS];
static uint32_t                         BfmSPI_BAT_id[NBRS];
static MUSPI_InjFifoSubGroup_t          BfmSPI_fifo_sub_group[NBRS];

static uint64_t            BfmSPI_HaloBytes[NBRS];


// Injection FIFO's (DMA senders)
static uint32_t                   BfmSPI_fifoids[NBRS] = {0,1,2,3,4,5,6,7};
static Kernel_InjFifoAttributes_t BfmSPI_injFifoAttrs[NBRS];


// Device memory resident buffers & reception counters
static void               *BfmSPI_send_mem;
static void               *BfmSPI_receive_mem;
static void               *BfmSPI_send_buffer[NBRS];
static void               *BfmSPI_receive_buffer[NBRS];
static void               *BfmSPI_simd_receive_buffer[2];
static volatile uint64_t  *BfmSPI_receive_counters;
// Put descriptors
static MUHWI_Descriptor_t *BfmSPI_descrs;
// What are these inj memory fifos????
static unsigned char      *BfmSPI_inj_memory_fifos[NBRS];

// Offsets from BasePA
static uint64_t  BfmSPI_offst_send_buffer[NBRS];
static uint64_t  BfmSPI_offst_receive_buffer[NBRS];
static uint64_t  BfmSPI_offst_receive_counters[NBRS];
static uint64_t  BfmSPI_offst_descr[NBRS];
static uint64_t  BfmSPI_offst_inj_memory_fifos[NBRS];
static uint      BfmSPI_HardwareThread;



// Little Dirac Op descriptors
static MUHWI_Descriptor_t *LdopSPI_descrs;
static volatile uint64_t  *LdopSPI_receive_counters;
static uint64_t  LdopSPI_offst_receive_counters[NBALL];
static uint64_t  LdopSPI_offst_descr[NBALL];
static uint64_t  LdopSPI_receive_buffer[NBALL];
static uint64_t  LdopSPI_send_buffer[NBALL];
static uint64_t  LdopSPI_offst_receive_buffer[NBALL];
static uint64_t  LdopSPI_offst_send_buffer[NBALL];

std::vector<MUHWI_Destination_t> LdopSPI_Destination;

void *comm_spi_malloc(size_t bytes)
{
  
  bytes = (bytes +64)&(~63UL);

  void * old = (void *)BfmSPI_device_malloc_ptr;
  BfmSPI_device_bytes     +=bytes;
  if ( BfmSPI_device_bytes > BFM_DEVICE_MEM_SIZE ) {
    printf("BAGEL: Out of SPI device memory\n");
    fflush(stdout);
    exit(-1);
  }
  BfmSPI_device_malloc_ptr+=bytes;
  return old;
}
void comm_spi_heap_restore(void)
{
  BfmSPI_device_malloc_ptr =BfmSPI_device_malloc_ptr_save;
  BfmSPI_device_bytes      =BfmSPI_device_bytes_save;
}
void comm_spi_heap_save(void)
{
  BfmSPI_device_malloc_ptr_save =BfmSPI_device_malloc_ptr;
  BfmSPI_device_bytes_save      =BfmSPI_device_bytes;
}

static int comm_spi_reserve = comm_spi_reserve_hardware(8);

int comm_spi_reserve_ldop_hardware(void);
int comm_spi_reserve_ldop_hardware(void)
{
  int rc;

  LdopSPI_descrs           = (MUHWI_Descriptor_t *) comm_spi_malloc(sizeof(MUHWI_Descriptor)*NBALL);
  LdopSPI_receive_counters = (uint64_t *)           comm_spi_malloc(sizeof(uint64_t)*NBALL);

  for(int ball=0;ball<NBALL;ball++){

    int dir = ball%NBRS;
    LdopSPI_offst_descr[ball]           = (uint64_t)&LdopSPI_descrs[ball]        -BfmSPI_BaseVa;
    LdopSPI_offst_receive_counters[ball]= MUSPI_GetAtomicOffsetFromBaseAddress (
									    &BfmSPI_BAT_sub_group[dir],
									     BfmSPI_BAT_id[dir],
									     (uint64_t)&LdopSPI_receive_counters[ball]-BfmSPI_BaseVa+BfmSPI_BasePa,
									     MUHWI_ATOMIC_OPCODE_STORE_ADD );

  }

  if(isBoss) {
    printf("bfmcommspi: Done ldop reservation\n");
    fflush(stdout);
  }

  return 0;
}

int comm_spi_reserve_hardware(int log2mb) 
{
  if ( BfmSPI_MU_resources_reserved ) return 1;

  Kernel_GetPersonality(&personality, sizeof(personality));
  isBoss = 1;
  if ( personality.Network_Config.Acoord!=0 ) isBoss=0;
  if ( personality.Network_Config.Bcoord!=0 ) isBoss=0;
  if ( personality.Network_Config.Ccoord!=0 ) isBoss=0;
  if ( personality.Network_Config.Dcoord!=0 ) isBoss=0;
  if ( personality.Network_Config.Ecoord!=0 ) isBoss=0;

  uint64_t rc;
  BFM_DEVICE_MEM_SIZE = (32*1024*1024) ;
  if ( log2mb == 6 ) BFM_DEVICE_MEM_SIZE = ( 64*1024*1024) ;
  if ( log2mb == 7 ) BFM_DEVICE_MEM_SIZE = (128*1024*1024) ;
  if ( log2mb == 8 ) BFM_DEVICE_MEM_SIZE = (256*1024*1024) ;
  if ( log2mb == 9 ) BFM_DEVICE_MEM_SIZE = (512*1024*1024) ;
  if ( log2mb ==10 ) BFM_DEVICE_MEM_SIZE = (1024*1024*1024) ;
  if ( log2mb > 10 ) if(isBoss) printf("Warning truncating the device memory requested to 1GB so don't be greedy\n");

  if(isBoss) printf("bfmcommspi: Printing all sorts of ugly stuff to convince you this is hard\n");
  if(isBoss) printf("bfmcommspi: SPI device memory pool %ld bytes\n",BFM_DEVICE_MEM_SIZE);

  ////////////////////////////////////////////////
  // Device memory 
  ////////////////////////////////////////////////
  BfmSPI_HardwareThread    = Kernel_PhysicalHWThreadID();
  BfmSPI_device_memory     = (unsigned char *)memalign(BFM_DEVICE_MEM_SIZE,BFM_DEVICE_MEM_SIZE);
  BfmSPI_device_malloc_ptr = BfmSPI_device_memory;
  BfmSPI_device_bytes      = 0;
  
  if(isBoss) printf("bfmcommspi: allocated SPI device memory pool in virtual memory space 0x%16.16lx - 0x%16.16lx \n",BfmSPI_device_memory,
			       (integer)BfmSPI_device_memory+BFM_DEVICE_MEM_SIZE);

  rc = Kernel_CreateMemoryRegion (&BfmSPI_device_memory_region,
				   BfmSPI_device_memory,
				   BFM_DEVICE_MEM_SIZE);
  CHECK("Kernel_CreateMemoryRegion\n");
  if(isBoss) printf("bfmcommspi: called Kernel_CreateMemoryRegion to pin this in contiguous physical region\n");

  BfmSPI_BasePa = (uint64_t)BfmSPI_device_memory_region.BasePa;
  BfmSPI_BaseVa = (uint64_t)BfmSPI_device_memory_region.BaseVa;

  if(isBoss) printf("bfmcommspi: Memory translation is as follows:\n");
  if(isBoss) printf("bfmcommspi: SPI device memory Base Physical Address %16.16lx\n",BfmSPI_BasePa);
  if(isBoss) printf("bfmcommspi: SPI device memory Base Virtual  Address %16.16lx\n",BfmSPI_BaseVa);
  
  for(int g=0;g<NBRS;g++){

    BfmSPI_BAT_id[g] = 0;

    rc = Kernel_AllocateBaseAddressTable (BfmSPI_sub_group[g],    // Thread/group
					 &BfmSPI_BAT_sub_group[g],// returned data struct
					  1,                   // 1 BAT id in list
					 &BfmSPI_BAT_id[g],    // Allocated entry
					  0 );                 // User use
    CHECK("Kernel_AllocateBaseAddressTable\n");

    if(isBoss) printf("bfmcommspi: allocating Block Address Translation (BAT) id %d Group %d\n",BfmSPI_BAT_id[g],BfmSPI_sub_group[g]);

    rc = MUSPI_SetBaseAddress (&BfmSPI_BAT_sub_group[g],
			        BfmSPI_BAT_id[g],
			        BfmSPI_BasePa );
    CHECK("MUSPI_SetBaseAddress\n");
  }

  
  /////////////////////////////////////////////////////////////////
  // Place Receive Counters, descriptors and buffers in device mem
  /////////////////////////////////////////////////////////////////
  if(isBoss) printf("bfmcommspi: preallocating reception counters, buffers, descriptors in this device memory pool\n");

  BfmSPI_descrs           = (MUHWI_Descriptor_t *) comm_spi_malloc(sizeof(MUHWI_Descriptor)*NBRS);
  BfmSPI_gsum_descr       = (MUHWI_Descriptor_t *) comm_spi_malloc(sizeof(MUHWI_Descriptor));
  BfmSPI_receive_counters = (uint64_t *)           comm_spi_malloc(sizeof(uint64_t)*NBRS);

  //////////////////////////////////////////////////////////////////////
  // Injection FIFO's
  //////////////////////////////////////////////////////////////////////
  if(isBoss)   printf("bfmcommspi: Injection FIFO setup\n");

  for(int i=0;i<NBRS;i++){
    BfmSPI_injFifoAttrs[i].RemoteGet = 0;
    BfmSPI_injFifoAttrs[i].System    = 0;
  }

  for(int i=0;i<NBRS;i++){
    uint32_t numfifo = 1;
    rc = Kernel_AllocateInjFifos (BfmSPI_sub_group[i], 
				  &BfmSPI_fifo_sub_group[i],
				  numfifo, 
				  &BfmSPI_fifoids[i],
				  &BfmSPI_injFifoAttrs[i]);

    CHECK("Kernel_AllocateInjFifos\n");
    
  }

  for(int dir=0;dir<NBRS;dir++){

    BfmSPI_inj_memory_fifos[dir]= (unsigned char *)comm_spi_malloc(INJ_MEMORY_FIFO_SIZE);
    
    BfmSPI_offst_inj_memory_fifos[dir]= (uint64_t)BfmSPI_inj_memory_fifos[dir]-BfmSPI_BaseVa;
    BfmSPI_offst_descr[dir]           = (uint64_t)&BfmSPI_descrs[dir]         -BfmSPI_BaseVa;
    BfmSPI_offst_receive_counters[dir] = MUSPI_GetAtomicOffsetFromBaseAddress ( 
									    &BfmSPI_BAT_sub_group[dir],
									     BfmSPI_BAT_id[dir],
									     (uint64_t)&BfmSPI_receive_counters[dir]
									    -BfmSPI_BaseVa
									    +BfmSPI_BasePa,
									     MUHWI_ATOMIC_OPCODE_STORE_ADD );


    rc = Kernel_InjFifoInit (&BfmSPI_fifo_sub_group[dir],
			      BfmSPI_fifoids[dir],
			     &BfmSPI_device_memory_region,      
			      BfmSPI_offst_inj_memory_fifos[dir],
			      INJ_MEMORY_FIFO_SIZE-1);    
    CHECK("Kernel_InjFifoInit\n");
                
    rc = Kernel_InjFifoActivate (&BfmSPI_fifo_sub_group[dir], 1, &BfmSPI_fifoids[dir],
				 KERNEL_INJ_FIFO_ACTIVATE);    
    CHECK("Kernel_InjFifoActivate \n");

  }

  ///////////////////////////////////////////////////////////////
  // Fast barrier support -- MPI_Barrier appears to be "detuned"
  ///////////////////////////////////////////////////////////////
  rc = MUSPI_GIBarrierInit ( &BfmSPI_GIBarrier,BfmSPI_GIBarrierClassRoute);
  CHECK("MUSPI_GIBarrierInit\n");

  if(isBoss)   printf("bfmcommspi: Global Interrupt Barriers initialised\n");

  MUSPI_GIBarrierEnterAndWait ( &BfmSPI_GIBarrier );


  if(isBoss) printf("bfmcommspi: SPI setup is complete\n");
  if(isBoss) printf("bfmcommspi: End of scary messages\n");
  fflush(stdout);

  MUSPI_GIBarrierEnterAndWait ( &BfmSPI_GIBarrier );

  comm_spi_reserve_ldop_hardware();

  MUSPI_GIBarrierEnterAndWait ( &BfmSPI_GIBarrier );

  if(isBoss) printf("bfmcommspi: preallocated section of device heap state for reinitialisation\n");
  comm_spi_heap_save();

  return 0;
}

template <class Float> 
void bfmcommspi<Float>::comm_spi_barrier(void) 
{
#if 0
  int rc = MUSPI_GIBarrierEnterAndWait ( &BfmSPI_GIBarrier );
  CHECK("MUSPI_GIBarrierEnterAndWait\n");
  int smear = QMP_me.Destination.E_Destination;
  //  Delay(smear*2000);
#else 
  QMP_barrier();
#endif
}

template <class Float> void bfmcommspi<Float>::comm_gsum(double *vals,int N)
{
#if 0
  for(int i=0;i<N;i++){
    this->comm_gsum(vals[i]);
  }
#else
  int me = this->thread_barrier();
  int rc;

  if ( me == 0 ) { 

    // Preload reception counter for 1 word gsum
    BfmSPI_receive_counters[BfmSPI_gsum_buf] = N*sizeof(double);
    BfmSPI_gsum_descr->Message_Length=N*sizeof(double);
    for(int i=0;i<N;i++){
      ((double *)BfmSPI_receive_buffer[BfmSPI_gsum_buf])[i] = 0.0;    
      ((double *)BfmSPI_send_buffer[BfmSPI_gsum_buf])[i]    = vals[i];
    }
    ppc_msync();  // make sure we pick up the correct start for this round
    
    // Might not be necessary?
    comm_spi_barrier();

    // Send with SPI
    rc = MUSPI_InjFifoInject(MUSPI_IdToInjFifo( BfmSPI_fifoids[BfmSPI_gsum_buf],
					       &BfmSPI_fifo_sub_group[BfmSPI_gsum_buf]),
			                        BfmSPI_gsum_descr);
    if( rc == -1 ) { 
      this->Error("Gsum SPI FifoInject failed with code %d\n",rc);
      exit(-1);
    }

    // Wait till done
    uint32_t polls = 0;
    while( BfmSPI_receive_counters[BfmSPI_gsum_buf]!= 0 ) { 
      polls++;
      if ( polls == 1024*1024 ) {
	this->Error("Inordinate spinning in comm_gsum[vec] %d/%d\n",BfmSPI_receive_counters[BfmSPI_gsum_buf],N*sizeof(double));
	polls=0;
      }
    }
    ppc_msync();  // make sure we pick up the correct start for this round

    for(int i=0;i<N;i++){
      vals[i]=((double *)BfmSPI_receive_buffer[BfmSPI_gsum_buf])[i];
    }
    ppc_msync();
  }
  this->thread_barrier();
#endif
}

static int comm_spi_gsum_counter;

template <class Float> 
void bfmcommspi<Float>::comm_gsum(double &val)
{
  int me = this->thread_barrier();
  int rc;

  if ( me == 0 ) { 

    // Preload reception counter for 1 word gsum
    BfmSPI_receive_counters[BfmSPI_gsum_buf] = sizeof(double);
    BfmSPI_gsum_descr->Message_Length=sizeof(double);
    *((double *)BfmSPI_receive_buffer[BfmSPI_gsum_buf]) = 0.0;    
    *((double *)BfmSPI_send_buffer[BfmSPI_gsum_buf])    = val;

    ppc_msync();  // make sure we pick up the correct start for this round
    
    // Might not be necessary?
    //    comm_spi_barrier();

    // Send with SPI
    rc = MUSPI_InjFifoInject(MUSPI_IdToInjFifo( BfmSPI_fifoids[BfmSPI_gsum_buf],
						&BfmSPI_fifo_sub_group[BfmSPI_gsum_buf]),
			     BfmSPI_gsum_descr);
    if( rc == -1 ) { 
      this->Error("Gsum SPI FifoInject failed with code %d\n",rc);	fflush(stdout);
      exit(-1);
    }

    // Wait till done
    uint64_t polls = 0;
    int tid = omp_get_thread_num();
    comm_spi_gsum_counter++;
    while( BfmSPI_receive_counters[BfmSPI_gsum_buf]!= 0 ) { 
      polls++;
      if ( polls%( 1024*1024*1024) == (1024*1024*1024-1) ) {
	this->Error("Inordinate spinning tid-%d in gsum [%d]\n",tid,comm_spi_gsum_counter);
      }
    }
    ppc_msync();  // make sure we pick up the correct start for this round
  }

  this->thread_barrier();
  val = *((double *)BfmSPI_receive_buffer[BfmSPI_gsum_buf]);
  
}


template <class Float> 
void bfmcommspi<Float>::comm_init(void) 
{

  const int A=0;
  const int B=1;
  const int C=2;
  const int D=3;
  const int E=4;

  uint64_t t0,t1;
  
  /***************************
   *Reset the memory heap
   * - Warning coexistence of two bfms precludes single/double
   * - scope lock needed?
   ***************************
   */
  t0=GetTimeBase();
  comm_spi_heap_restore();
  t1=GetTimeBase();
  //  this->BossDebug("bfmcommspi: %s time %ld\n","heap_restore",t1-t0);

  /*
   * Allocate device memory
   */
  uint64_t words = this->simd_allbound*SPINOR_SIZE*this->simd()*this->cbLs*2;
  uint64_t bytes = 2*words * sizeof(Float);
  uint64_t mem_pool = (bytes + 1024*1024);

  if ( mem_pool > BFM_DEVICE_MEM_SIZE ) { 
    this->Error("SPI Memory pool exceeds reserved memory pool\n");exit(0);
  }

  //  t0=GetTimeBase();
  BfmSPI_send_mem    = comm_spi_malloc(words * sizeof(Float));
  BfmSPI_receive_mem = comm_spi_malloc(words * sizeof(Float));
  //  t1=GetTimeBase();
  //this->BossDebug("bfmcommspi: %s time %ld\n","malloc",t1-t0);

  t0=GetTimeBase();
  for(int mu=0;mu<4;mu++){
    for(int pm=0;pm<2;pm++){

      int dir= pm+2*mu;
      int simd_factor = this->simd();
      //      simd_factor = simd_factor/this->simd_grid[mu];
      int words = this->simd_nbound[mu] * HALF_SPINOR_SIZE * simd_factor * this->cbLs;

      BfmSPI_send_buffer[dir]    = &( (Float *)BfmSPI_send_mem)  [this->comm_offset[pm][mu]*SPINOR_SIZE*this->simd()*this->cbLs];
      BfmSPI_receive_buffer[dir] = &( (Float *)BfmSPI_receive_mem)[this->comm_offset[pm][mu]*SPINOR_SIZE*this->simd()*this->cbLs];

      if (this->simd_grid[mu] > 1 ){
	uint64_t simd_words = this->simd_nbound[mu] * SPINOR_SIZE * this->cbLs * 2;
	uint64_t simd_bytes = words * sizeof(Float);
	BfmSPI_simd_receive_buffer[pm]    = comm_spi_malloc(simd_bytes);
	BfmSPI_offst_receive_buffer[dir]  = (uint64_t)BfmSPI_simd_receive_buffer[pm]  -BfmSPI_BaseVa;
      } else {
	BfmSPI_offst_receive_buffer[dir]  = (uint64_t)BfmSPI_receive_buffer[dir]  -BfmSPI_BaseVa;
      }

      // Virt->Phys
      BfmSPI_offst_send_buffer[dir]     = (uint64_t)BfmSPI_send_buffer[dir]     -BfmSPI_BaseVa;

    }
  }
  t1=GetTimeBase();
  //this->BossDebug("bfmcommspi: %s time %ld\n","offsets",t1-t0);

  t0=GetTimeBase();
  for(int mu=0;mu<4;mu++){
    for(int pm=0;pm<2;pm++){

      int dir= pm+2*mu;
      int simd_factor = this->simd();
      simd_factor = simd_factor/this->simd_grid[mu];
      int words = this->simd_nbound[mu] * HALF_SPINOR_SIZE * simd_factor * this->cbLs;

      BfmSPI_HaloBytes[dir]       = words*sizeof(Float);      
      //    this->BossDebug("bfmcommspi: Halo bytes[%d] = %d\n",dir,BfmSPI_HaloBytes[dir]);
      this->sendbufs[dir] = (Float *)BfmSPI_send_buffer[dir];
      if ( (this->simd_grid[mu] > 1) && !this->local_comm[mu]) {
	this->recvbufs[dir]  = (Float *)BfmSPI_receive_buffer[dir];
	this->simd_rbuf[dir] = (Float *)BfmSPI_simd_receive_buffer[pm];
      } else {
	this->recvbufs[dir] = (Float *)BfmSPI_receive_buffer[dir];
      }

    }
  }  
  t1=GetTimeBase();
  //  this->BossDebug("bfmcommspi: %s time %ld\n","buffers",t1-t0);

  ////////////////////////////////////////////////
  // DMA/FIFO controls
  ////////////////////////////////////////////////

  {
  MUSPI_Pt2PtDirectPutDescriptorInfo_t info;
  MUSPI_Pt2PtDescriptorInfoFields_t Pt2Pt;/**< Point-to-point descriptor info */
  MUSPI_DirectPutDescriptorInfoFields_t DirectPut;

  
  /*
  this->BossDebug("bfmcommspi: DMA/FIFO controls\n");
  this->BossDebug("bfmcommspi: Using Pt2PtDirectPut messages\n");
  this->BossDebug("bfmcommspi: Using SPIGIBarrier() to synchronize\n");
  this->BossDebug("bfmcommspi: Retrieving neigbour information from QMP and checking it matches hardware torus\n");
  this->BossDebug("bfmcommspi: Must use PAB's QMP hack -qmp-geom native/native3d \n");
  */

  t0=GetTimeBase();
  for(int mu=0;mu<4;mu++){
    for(int pm=0;pm<2;pm++){

      int dir = pm+2*mu;
      int buf= (1-pm)+2*mu;

      int mach_dim = -1;
      int ndelta   = 0;
      int L=0;
      int me=0;
      int you=0;

      if ( QMP_neigh[dir].Destination.Destination != QMP_me.Destination.Destination ) {

	if ( QMP_neigh[dir].Destination.A_Destination != QMP_me.Destination.A_Destination ) {
	  mach_dim=A; ndelta++;
	}
	if ( QMP_neigh[dir].Destination.B_Destination != QMP_me.Destination.B_Destination ) {
	  mach_dim=B; ndelta++;
	}
	if ( QMP_neigh[dir].Destination.C_Destination != QMP_me.Destination.C_Destination ) {
	  mach_dim=C; ndelta++;
	}
	if ( QMP_neigh[dir].Destination.D_Destination != QMP_me.Destination.D_Destination ) {
	  mach_dim=D; ndelta++;
	}
	if ( QMP_neigh[dir].Destination.E_Destination != QMP_me.Destination.E_Destination ) {
	  mach_dim=E; ndelta++;
	}
	if ( ndelta!= 1) { 
	  this->Error("bfmspi: mapping logic bomb 2\n");
	}

	switch(mach_dim){ 
	case A:
	  me =QMP_me.Destination.A_Destination;
	  you=QMP_neigh[dir].Destination.A_Destination;
	  L  =QMP_Labcde[A];
	  break;
	case B:
	  me =QMP_me.Destination.B_Destination;
	  you=QMP_neigh[dir].Destination.B_Destination;
	  L  =QMP_Labcde[B];
	  break;
	case C:
	  me =QMP_me.Destination.C_Destination;
	  you=QMP_neigh[dir].Destination.C_Destination;
	  L  =QMP_Labcde[C];
	  break;
	case D:
	  me =QMP_me.Destination.D_Destination;
	  you=QMP_neigh[dir].Destination.D_Destination;
	  L  =QMP_Labcde[D];
	  break;
	case E:
	  me =QMP_me.Destination.E_Destination;
	  you=QMP_neigh[dir].Destination.E_Destination;
	  L  =QMP_Labcde[E];
	  break;
	default:
	  this->Error("bfmspi: mapping logic bomb 3\n");
	  exit(-1);
	}
	
	int mach_pm=-1;
	if ( me == ((you+1)) ) { 
	  mach_pm=0;
	} else if ( you = ((me+1)) ){
	  mach_pm=1;
	} else if ( me == ((you+1)%L) ) { // Avoid peri wrap if possible
	  mach_pm=0;
	} else if ( you = ((me+1)%L) ){
	  mach_pm=1;
	} else { 
	  this->Error("bfmspi: mapping logic bomb 4\n");
	  exit(-1);
	}
	  
 	int mach_dir = mach_dim*2 + mach_pm;

        uint64_t flagbit = FW_BIT(mach_dir);
	if ( !(personality.Network_Config.NetFlags2 & flagbit) ) {
	  this->Error("machdir %d is not enabled\n",mach_dir);
	}
	
	info.Base.Pre_Fetch_Only = 0;
	info.Base.Payload_Address= BfmSPI_offst_send_buffer[buf]+BfmSPI_BasePa; // phys addr
	info.Base.Message_Length = BfmSPI_HaloBytes[buf];  
	info.Base.Torus_FIFO_Map = BfmSPI_Torus_FIFO_Map[mach_dir];
	info.Base.Dest           = QMP_neigh[dir];

	///////////////////////////
	// Debug hack cache values
	///////////////////////////
	last_logical_dir=dir;

	info.Pt2Pt.Hints_ABCD = BfmSPI_hints_abcd[mach_dir];
	info.Pt2Pt.Misc1 = BfmSPI_hints_e[mach_dir];
	info.Pt2Pt.Misc2 = MUHWI_PACKET_VIRTUAL_CHANNEL_DETERMINISTIC;
	info.Pt2Pt.Skip  = 8;

	info.DirectPut.Rec_Payload_Base_Address_Id = BfmSPI_BAT_id[dir];
	info.DirectPut.Rec_Payload_Offset          = BfmSPI_offst_receive_buffer[buf];
	info.DirectPut.Rec_Counter_Base_Address_Id = BfmSPI_BAT_id[dir];
	info.DirectPut.Rec_Counter_Offset          = BfmSPI_offst_receive_counters[0];
	info.DirectPut.Pacing                      = MUHWI_PACKET_DIRECT_PUT_IS_NOT_PACED;

 	int rc = MUSPI_CreatePt2PtDirectPutDescriptor(&BfmSPI_descrs[dir],&info);
	CHECK("MUSPI_CreatePt2PtDirectPutDescriptor \n");
	
	//	if ( this->isBoss() ) {
	//	  MUSPI_DescriptorDumpHex("Dslash packet...",&BfmSPI_descrs[dir]);fflush(stdout);
	//	}

      } else {
	if ( !this->local_comm[mu] ) {
	  this->Error("bfmspi: mapping logic bomb 1\n");
	  exit(-1);
	}
      }
    }
  }
  }
  t1=GetTimeBase();

  /*
  this->BossDebug("bfmcommspi: %s time %ld\n","DirectPut",t1-t0);
  this->BossDebug("bfmcommspi: Created all SPI Pt2PtDirectPut descriptors \n");
  */
  ///////////////////////////////////////////
  // Fast global sum
  ///////////////////////////////////////////
  t0=GetTimeBase();
  {
    
    MUHWI_Destination_t dest;
    MUSPI_SetUpDestination ( &dest, 0, 0, 0, 0, 0 );

    MUSPI_CollectiveDirectPutDescriptorInfo_t  info;

    // Reuse send buffer [0] and FIFO_Map[0] for the send
    info.Base.Pre_Fetch_Only = 0;
    info.Base.Payload_Address= BfmSPI_offst_send_buffer[BfmSPI_gsum_buf]+BfmSPI_BasePa; // phys addr
    info.Base.Message_Length = sizeof(double);
    info.Base.Torus_FIFO_Map = MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_CUSER;
    info.Base.Dest           = dest;

    info.DirectPut.Rec_Payload_Base_Address_Id = BfmSPI_BAT_id[BfmSPI_gsum_buf];
    info.DirectPut.Rec_Payload_Offset          = BfmSPI_offst_receive_buffer[BfmSPI_gsum_buf];
    info.DirectPut.Rec_Counter_Base_Address_Id = BfmSPI_BAT_id[BfmSPI_gsum_buf];
    info.DirectPut.Rec_Counter_Offset          = BfmSPI_offst_receive_counters[BfmSPI_gsum_buf];
    info.DirectPut.Pacing                      = MUHWI_PACKET_DIRECT_PUT_IS_NOT_PACED;
  
    info.Collective.Class_Route = 0;
    info.Collective.Word_Length = sizeof(double);
    info.Collective.Op_Code     = MUHWI_COLLECTIVE_OP_CODE_FLOATING_POINT_ADD;
    info.Collective.Misc        = MUHWI_PACKET_VIRTUAL_CHANNEL_USER_SUB_COMM | MUHWI_COLLECTIVE_TYPE_ALLREDUCE;
    info.Collective.Skip        = 0;

    int rc = MUSPI_CreateCollectiveDirectPutDescriptor(BfmSPI_gsum_descr,&info);
    CHECK("MUSPI_CreateCollectiveDirectPutDescriptor \n");



  }

  allBytes=0;
  for(int mu=0;mu<4;mu++){
    if ( !this->local_comm[mu] ) {
      allBytes+=2*BfmSPI_HaloBytes[2*mu];  
    }
  }
  t1=GetTimeBase();

  haloes.resize(0);
  haloes.push_back(make_pair<int,int>(BfmSPI_HaloBytes[0],0));  
  haloes.push_back(make_pair<int,int>(BfmSPI_HaloBytes[2],1));  
  haloes.push_back(make_pair<int,int>(BfmSPI_HaloBytes[4],2));  
  haloes.push_back(make_pair<int,int>(BfmSPI_HaloBytes[6],3));  
  std::sort(haloes.rbegin(),haloes.rend());

  /*
  this->BossDebug("bfmcommspi: %s time %ld\n","Collective",t1-t0);
  */  
  QMP_barrier();
  comm_spi_barrier();
}
template <class Float> 
void bfmcommspi<Float>::comm_end(void) 
{
  // We simply reinit everything in static allocation
}

template <class Float> 
void bfmcommspi<Float>::comm_start_precision_test(void)
{
  if ( sizeof(Float)==4 ) {
    for(int mu=3;mu<4;mu++){
      for(int pm=0;pm<2;pm++){
	int dir= pm+2*mu;
	int simd_factor = this->simd();
	simd_factor = simd_factor/this->simd_grid[mu];

	int words = this->simd_nbound[mu] * HALF_SPINOR_SIZE * simd_factor * this->cbLs; //N2SPIN*NCOLOR*NREIM
	uint16_t *buf = (uint16_t *)this->sendbufs[dir];
	uint32_t conv;
	for(int w=0;w<32;w+=8){
	//	for(int w=0;w<words*2;w+=8){
	  ThreadBossLog("%d %d %d : %4.4x %4.4x %4.4x %4.4x : %4.4x %4.4x %4.4x %4.4x\n",
			mu,pm,w,
                        buf[w+0], buf[w+1], buf[w+2], buf[w+3],
                        buf[w+4], buf[w+5], buf[w+6], buf[w+7]);

	}
      }
    }
  }
}

template <class Float> 
void bfmcommspi<Float>::comm_start(int result_cb, Fermion_t psi,int dag) 
{
  // gather the faces. Routines here are threaded.
  // All threads cooperate in gathering.
  static int printed;

  int me = this->thread_barrier();

  if ( !me ) {
    int len = allBytes;
    if ( this->precision_test ) len = len/2;
    BfmSPI_receive_counters[0]=len;
    ppc_msync();
    this->comm_spi_barrier();
  }
  //  this->thread_barrier();


  uint64_t t_gat =0;
  uint64_t t_spi =0;
  uint64_t t0=0;

  for(int mup=0;mup<4;mup++){

    int mu = haloes[mup].second;

    t0=GetTimeBase();
    this->gather (mu,result_cb,psi,dag);
    t_gat+=GetTimeBase()-t0;

    this->thread_barrier();

    t0=GetTimeBase();
    comm_spi_start(mu);
    t_spi+=GetTimeBase()-t0;
  }

  if ( me == 0 ) {
    //    if(this->iter == this->time_report_iter+2){ 
    if ( this->precision_test ){
      if (printed==9) { 
	this->BossLog("gather    : %d cyc\n",t_gat);
	this->BossLog("spi_start : %d cyc\n",t_spi);
      }    
      if ( this->precision_test )  printed++;
    }
  }
  this->thread_barrier();
  
  return;
}

template <class Float> 
void bfmcommspi<Float>::comm_complete(int result_cb, Fermion_t psi) 
{
  static int printed;
  int me=this->thread_barrier();
  uint64_t t1=GetTimeBase();
  if ( me == 0 ) {
    comm_spi_complete() ;
  }
  this->thread_barrier();
  uint64_t t2=GetTimeBase();
  this->comm_merge();
  this->thread_barrier();
  uint64_t t3=GetTimeBase();
  if ( me == 0 ) {
    if ( this->precision_test ){
      if (printed==5) { 
	this->BossLog("spi_complete    : %d cyc\n",t2-t1);
	this->BossLog("merge           : %d cyc\n",t3-t2);
      }    
      if ( this->precision_test )  printed++;
    }
  }
  return;
}
template <class Float> 
void bfmcommspi<Float>::comm_spi_start(int mu) 
{
  int me, throff,thrlen;
  int npm=2;
  this->thread_work(npm,me,thrlen,throff);


  for(int pm=throff;pm<throff+thrlen;pm++){

    int dir = pm + 2*mu;
    if ( !this->local_comm[mu] ) {

      int len = BfmSPI_HaloBytes[dir];  

      if ( this->precision_test ) len = len/2;
      BfmSPI_descrs[dir].Message_Length=len;
      ppc_msync();

      uint64_t rc = MUSPI_InjFifoInject(MUSPI_IdToInjFifo( BfmSPI_fifoids[dir],&BfmSPI_fifo_sub_group[dir]),&BfmSPI_descrs[dir]);
      if ( rc == ((uint64_t)-1) ) {
	this->Error("MUSPI_InjFifoInject\n");
	exit(-1);
      }
    }
  }
  this->thread_barrier();
}

template <class Float> 
void bfmcommspi<Float>::comm_spi_complete(void) 
{
  uint64_t old = 0;
  uint64_t stale=0;
  uint64_t polls=0;

  while( BfmSPI_receive_counters[0]!= 0 ) { 
    polls++;
    if ( polls == 1024*1024*1024 ) {
      polls=0;
      this->Error("Inordinate spinning in comms receive\n");fflush(stdout);
    }
  }
  ppc_msync();
  this->comm_spi_barrier();
}

///////////////////////////////////////
// Little dirac operator support
///////////////////////////////////////
#include <BfmHDCG.h>


template<class cFloat>
void HaloInitBuffers(BfmHDCG<cFloat> *BMG)
{
  int Nmsg = BMG->tproc.size();

  int mynode = QMP_get_node_number();
  int  nodes = QMP_get_number_of_nodes();

  // Build map of QMP node # to MU destination 
  LdopSPI_Destination.resize(nodes);
  for(int i=0;i<nodes;i++){
    int dest=0;
    if ( i==mynode ) {
      dest=QMP_me.Destination.Destination;
    }
    QMP_sum_int(&dest);
    LdopSPI_Destination[i].Destination.Destination=(uint32_t)dest;
  }

  ////////////////////////////////////////////////
  // Buffers
  ////////////////////////////////////////////////
  if ( !BfmHDCGStatic::StaticInitialised ) { 

    BfmHDCGStatic::sendbufs_static.resize(Nmsg);
    BfmHDCGStatic::recvbufs_static.resize(Nmsg);

    int pad=1;
    if ( sizeof(cFloat) == 4 ) pad = 2;

    for(int m=0;m<Nmsg;m++){
      BfmHDCGStatic::sendbufs_static[m] = (std::complex<double> *)comm_spi_malloc(BMG->sendbuf_bytes[m]*pad);
      LdopSPI_send_buffer[m]   =(uint64_t)BfmHDCGStatic::sendbufs_static[m];
    }
    
    for(int m=0;m<Nmsg;m++){
      BfmHDCGStatic::recvbufs_static[m] = (std::complex<double> *)comm_spi_malloc(BMG->recvbuf_bytes[m]*pad);
      LdopSPI_receive_buffer[m]=(uint64_t)BfmHDCGStatic::recvbufs_static[m];
    }

    for(int m=0;m<Nmsg;m++){
      // Virt->Phys
      LdopSPI_offst_receive_buffer[m] = (uint64_t)LdopSPI_receive_buffer[m]  -BfmSPI_BaseVa;
      LdopSPI_offst_send_buffer[m]    = (uint64_t)LdopSPI_send_buffer[m]     -BfmSPI_BaseVa;
    }
    BfmHDCGStatic::StaticInitialised=1;
  }

  for(int m=0;m<Nmsg;m++){
    BMG->sendbufs[m] = (std::complex<cFloat> *)BfmHDCGStatic::sendbufs_static[m];
    BMG->recvbufs[m] = (std::complex<cFloat> *)BfmHDCGStatic::recvbufs_static[m];
  }

  ////////////////////////////////////////////////
  // DMA/FIFO controls
  ////////////////////////////////////////////////
  BMG->linop_d->comm_spi_barrier();

  {
    MUSPI_Pt2PtDirectPutDescriptorInfo_t info;
    MUSPI_Pt2PtDescriptorInfoFields_t Pt2Pt;/**< Point-to-point descriptor info */
    MUSPI_DirectPutDescriptorInfoFields_t DirectPut;

    if (BMG->linop_d->isBoss()) {
      BMG->linop_d->BossDebug("BfmHDCGSPI: DMA/FIFO controls\n");
      BMG->linop_d->BossDebug("BfmHDCGSPI: Using Pt2PtDirectPut messages\n");
      fflush(stdout);
    }

    for(int m=0;m<Nmsg;m++){

      //      BMG->linop_d->BossDebug("BfmHDCGSPI: Message %d size %d to proc %d\n",
      // m,BMG->sendbuf_bytes[m],BMG->tproc[m]	       );

      BMG->linop_d->comm_spi_barrier();
	
      //      FIXME MUST PICK A MACH_DIR
	info.Base.Pre_Fetch_Only = 0;
	info.Base.Payload_Address= LdopSPI_offst_send_buffer[m]+BfmSPI_BasePa; // phys addr
	info.Base.Message_Length = BMG->sendbuf_bytes[m];
	info.Base.Dest           = LdopSPI_Destination[BMG->tproc[m]];

	// Any fifo
	info.Base.Torus_FIFO_Map = 
	    MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_AM
	  | MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_AP
	  | MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_BM
	  | MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_BP
	  | MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_CM
	  | MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_CP
	  | MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_DM
	  | MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_DP
	  | MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_EM
	  | MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_EP;

	// No hints
	info.Pt2Pt.Hints_ABCD = MUHWI_PACKET_HINT_A_NONE;
	info.Pt2Pt.Misc1      = MUHWI_PACKET_HINT_E_NONE|MUHWI_PACKET_DO_NOT_ROUTE_TO_IO_NODE|MUHWI_PACKET_USE_DETERMINISTIC_ROUTING| MUHWI_PACKET_DO_NOT_DEPOSIT;
	info.Pt2Pt.Misc2      = MUHWI_PACKET_VIRTUAL_CHANNEL_DETERMINISTIC;
	info.Pt2Pt.Skip       = 8;

	info.DirectPut.Rec_Payload_Base_Address_Id = BfmSPI_BAT_id[0];
	info.DirectPut.Rec_Counter_Base_Address_Id = BfmSPI_BAT_id[0];
	info.DirectPut.Rec_Payload_Offset          = LdopSPI_offst_receive_buffer[m];
	info.DirectPut.Rec_Counter_Offset          = LdopSPI_offst_receive_counters[0];

	info.DirectPut.Pacing                      = MUHWI_PACKET_DIRECT_PUT_IS_NOT_PACED;
	
 	int rc = MUSPI_CreatePt2PtDirectPutDescriptor(&LdopSPI_descrs[m],&info);
	CHECK("MUSPI_CreatePt2PtDirectPutDescriptor \n");

	//	if (BMG->linop->isBoss()) {
	//	  MUSPI_DescriptorDumpHex("HDCG packet...",&LdopSPI_descrs[m]);fflush(stdout);
	//	}

    }
  }  
}
template<class cFloat>
void HaloExchangeCommComplete(BfmHDCG<cFloat> *BMG, int depth)
{
  int me = BMG->linop_d->thread_barrier();
  uint64_t old = 0;
  uint64_t stale=0;
  uint64_t polls=0;
  if ( me == 0 ) { 
    while( LdopSPI_receive_counters[0]!= 0 ) { 
      polls++;
      if ( polls == 1024*1024*1024 ) {
	polls=0;
	BMG->linop_d->Error(
			    "Inordinate spinning in comms receive for little dirac op: counter 0x%lx\n",
			    (uint64_t)LdopSPI_receive_counters[0]);
	fflush(stdout);
      }
    }
    //    ppc_msync();
  }
  BMG->linop_d->thread_barrier();
}

template<class cFloat>
void HaloExchangeCommStart(BfmHDCG<cFloat> *BMG,int depth)
{
  int me,msgnum,msgoff; 
  int Nmsg = BMG->tproc.size();
  int incr=0;

  int throff,thrlen;

  BMG->linop_d->thread_work(NBRS,me,thrlen,throff);

  if(me==0) {
    for(int m=0;m<Nmsg;m++) {
      //      incr+=BMG->sendbuf_bytes[m];
      incr+=BMG->sendbuf_depth_bytes[depth][m];
    }

    LdopSPI_receive_counters[0]=incr;
    ppc_msync();
    BMG->linop_d->comm_spi_barrier();
  }

  BMG->linop_d->thread_barrier();
  
  // Put a different thread onto each of the 8 FIFO's allocated by BfmCommSpi
  // and stuff the packets in
  for(int dir=throff;dir<throff+thrlen;dir++){

    int num_msg,first_msg;
    ThreadModelSingle WorkAllocator; 
    WorkAllocator.nthread=NBRS;
    WorkAllocator.thread_work_nobarrier(Nmsg,dir,num_msg,first_msg);

#if 1
    for(int msg=first_msg;msg<first_msg+num_msg;msg++){
      //      LdopSPI_descrs[msg].Message_Length=BMG->sendbuf_bytes[msg];
      int len = BMG->sendbuf_depth_bytes[depth][msg];
      if ( len ) {
	LdopSPI_descrs[msg].Message_Length=len;
	uint64_t rc = MUSPI_InjFifoInject(MUSPI_IdToInjFifo( BfmSPI_fifoids[dir],&BfmSPI_fifo_sub_group[dir]),&LdopSPI_descrs[msg]);
	
	if ( rc == ((uint64_t)-1) ) {
	  BMG->linop_d->Error("MUSPI_InjFifoInject\n");fflush(stdout);
	  exit(-1);
	}
      //      while(!MUSPI_CheckDescComplete (MUSPI_IdToInjFifo( BfmSPI_fifoids[dir],&BfmSPI_fifo_sub_group[dir]),rc )) { };
      }
    }
#else
    for(int msg=first_msg;msg<first_msg+num_msg;msg++){
      LdopSPI_descrs[msg].Message_Length=BMG->sendbuf_bytes[msg];
    }
    uint64_t rc= MUSPI_InjFifoInjectMultiple (MUSPI_IdToInjFifo( BfmSPI_fifoids[dir],&BfmSPI_fifo_sub_group[dir]),&LdopSPI_descrs[first_msg],num_msg);

    if ( rc == ((uint64_t)-1) ) {
      BMG->linop_d->Error("MUSPI_InjFifoInject\n");fflush(stdout);
      exit(-1);
    }
#endif
  }

  BMG->linop_d->thread_barrier();
}

template void HaloInitBuffers(BfmHDCG<float> *BMG);
template void HaloInitBuffers(BfmHDCG<double> *BMG);
template void HaloExchangeCommComplete(BfmHDCG<float> *BMG,int depth);
template void HaloExchangeCommComplete(BfmHDCG<double> *BMG,int depth);
template void HaloExchangeCommStart(BfmHDCG<float> *BMG,int depth);
template void HaloExchangeCommStart(BfmHDCG<double> *BMG,int depth);

template class bfmcommspi <float>;
template class bfmcommspi <double>;
 
