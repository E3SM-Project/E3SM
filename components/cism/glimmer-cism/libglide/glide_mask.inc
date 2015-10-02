#ifndef GLIDE_MASK_INC

!Bits 1:0 - Ice presence (0 if no ice). NOTE: If bit 2 is activated, bit 1 must be activated
#define GLIDE_MASK_HAS_ICE             1
#define GLIDE_MASK_THIN_ICE            3
#define GLIDE_ICE_PRESENCE_BITS        3

!TODO - Is the grounding line mask needed?
!       Also, note that MASK_LAND is redundant given MASK_OCEAN, since one is the complement of the other.

!Bits 4:2 - Type of base (Land or ocean - grounding line has both bits on). The 16 bit specifies grounding line, 
! set up so that those points are treated as grounded)
#define GLIDE_MASK_LAND                4
#define GLIDE_MASK_OCEAN               8
#define GLIDE_MASK_GROUNDING_LINE      16 

!Bit 5: Identifies a margin (jump from zero to nonzero thickness). Margin type determined by whether ice grounded or floating.
#define GLIDE_MASK_MARGIN 32

!Bit 6: Identifies a dirichlet condition. The velocity at points marked with this bit should be held constant.
#define GLIDE_MASK_DIRICHLET_BC 64

!Bit 7: Identifies a computational domain boundary.  These are normally just activated on the edges of the domain,
!unless there is a domain decomposition (in which case they may be missing)
#define GLIDE_MASK_COMP_DOMAIN_BND 128 

!=======
! All mask values actually used in code (as defined below) should be made up of some combination of one or more 
! of the above "base" type bits.  

!TODO - Rename to GLIDE_ICEFREE_OCEAN?
!Checks for an iceless square
!Checks for open ocean with no ice
#define GLIDE_IS_OCEAN(mask)       (iand(mask, GLIDE_MASK_OCEAN) == GLIDE_MASK_OCEAN .and. GLIDE_NO_ICE(mask))

!TODO - Rename to GLIDE_ICEFREE_LAND?
!Checks for land with no ice
#define GLIDE_IS_LAND(mask)        (iand(mask, GLIDE_MASK_LAND) == GLIDE_MASK_LAND .and. GLIDE_NO_ICE(mask))

!Checks for the presence of any ice, dynamic or not
#define GLIDE_HAS_ICE(mask)        (iand(mask, GLIDE_MASK_HAS_ICE) == GLIDE_MASK_HAS_ICE)

!Checks for a lack of ice
#define GLIDE_NO_ICE(mask)         (iand(mask, GLIDE_MASK_HAS_ICE) == 0)

!Checks for the presence of ice that is below the ice dynamics limit
#define GLIDE_IS_THIN(mask)          (iand(mask,GLIDE_MASK_THIN_ICE)  == GLIDE_MASK_THIN_ICE)

!Checks for any ice, dynamic or not, that is on an ice shelf.
#define GLIDE_IS_FLOAT(mask)         (iand(mask,GLIDE_MASK_OCEAN) == GLIDE_MASK_OCEAN .and. GLIDE_HAS_ICE(mask))

!Checks for any ice, dynamic or not, that is grounded
#define GLIDE_IS_GROUND(mask)      (iand(mask,GLIDE_MASK_LAND) == GLIDE_MASK_LAND .and. GLIDE_HAS_ICE(mask))

!TODO - Remove hardwiring?
!       Actually, not sure this is needed
!Checks for any ice, dynamic or not, that is on the grounding line

!TODO: define from above available combinations rather than hardcode?
!17 = 16 + 1
#define GLIDE_IS_GROUNDING_LINE(mask) (iand(mask, 17) == 17)

!TODO - This one probably is not needed - NOTE that removing it affects the iarea/ivol calculation
!Checks for any ice, dynamic or not, that is either floating *or* on the grounding line
!TODO: IF NEEDED, define from above available combinations rather than hardcode?
#define GLIDE_IS_FLOAT_OR_GNDLINE(mask) (iand(mask, 24) > 0 .and. GLIDE_HAS_ICE(mask))

!TODO - This one probably is not needed - NOTE that removing it affects the iarea/ivol calculation
!Checks for any ice, dynamic or not, that is either grounded *or* on the grounding line
!TODO: IF NEEDED, define from above available combinations rather than hardcode?
#define GLIDE_IS_GROUND_OR_GNDLINE(mask) (iand(mask, 20) > 0 .and. GLIDE_HAS_ICE(mask))

!Checks whether this is an ice margin (thickness jumps from 0 to non-zero at this point)
#define GLIDE_IS_MARGIN(mask) (iand(mask, GLIDE_MASK_MARGIN) == GLIDE_MASK_MARGIN)

!TODO - Not sure this is needed
!Checks whether this is a margin in contact with the ocean, floating or not
#define GLIDE_IS_MARINE_ICE_EDGE(mask)    (GLIDE_IS_MARGIN(mask) .and. GLIDE_IS_FLOAT_OR_GNDLINE(mask))

!TODO - Not a good name for this mask
!Checks whether this is a margin in contact with the ocean
!41 = 32 + 8 + 1
!TODO: define from above available combinations rather than hardcode?
#define GLIDE_IS_CALVING(mask)            (iand(mask, 41) == 41)

!Checks whether this is a land margin
!37 = 32 + 4 + 1
!TODO: define from above available combinations rather than hardcode?
#define GLIDE_IS_LAND_MARGIN(mask)        (iand(mask, 37) == 37)

!TODO - Where are the Dirichlet and domain_bnd masks set in the code?
!Checks whether a dirichlet boundary has been defined at this point
#define GLIDE_IS_DIRICHLET_BOUNDARY(mask) (iand(mask, GLIDE_MASK_DIRICHLET_BC) == GLIDE_MASK_DIRICHLET_BC)

!Checks whether we are at the edge of the computational domain *and* there is ice in this square

!TODO: define from above available combinations rather than hardcode?

! 129 = 128 + 1
#define GLIDE_IS_COMP_DOMAIN_BND(mask)    (iand(mask, 129) == 129)

! table of common combinations:
!1			! has ice
!3			! has thin ice
!4		 	! land
!5 = 4 + 1              ! grounded ice
!7 = 4 + 3              ! grounded thin ice
!8		 	! ocean
!9 = 8 + 1              ! floating ice
!11 = 8 + 3             ! floating thin ice
!17 = 16 + 1		! grounding line with ice
!19 = 16 + 3		! grounding line with thin ice
!21 = 16 + 4 + 1	! non-thin ice grounding line
!23 = 16 + 4 + 3	! thin ice grounding line
!32			! margin
!37 = 32 + 4 + 1	! land margin
!39 = 32 + 4 + 3	! land margin with thin ice
!41 = 32 + 8 + 1	! calving margin (ocean)
!43 = 32 + 8 + 3	! calving margin (ocean) with thin ice
!53 = 32 + 21		! margin AND grounding line   
!55 = 32 + 23		! margin AND grounding line with thin ice  

#endif
