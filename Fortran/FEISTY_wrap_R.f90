! for R
subroutine f_setupbasic(pprod,bprod,T)
   use FEISTY !, only:
   use globals

   real(dp), intent(in)::pprod,bprod,T

   call setupbasic(pprod,bprod,T)
end subroutine f_setupbasic

subroutine f_setupbasic2(pprod, bprod, nStages, T)
   use FEISTY ! only : setupbasic2
   use globals

   real(dp), intent(in):: pprod, bprod, T
   integer, intent(in):: nStages

   call setupbasic2(pprod, bprod, nStages, T)
end subroutine f_setupbasic2

subroutine f_setupvertical(pprod, nStages,region)
   use FEISTY ! only :
   use globals

   real(dp), intent(in):: pprod
   integer, intent(in):: nStages,region

   call setupvertical(pprod, nStages,region)
end subroutine f_setupvertical

subroutine f_setupsquid(pprod, bottom, nStages)
   use FEISTY ! only :
   use globals

   real(dp), intent(in):: pprod, bottom
   integer, intent(in) :: nStages

   call setupsquid(pprod, bottom, nStages)
end subroutine f_setupsquid

subroutine f_calcderivatives(u, dudt)
   use FEISTY !, only:
   use globals

   real(dp), intent(in):: u(nGrid)
   real(dp), intent(out):: dudt(nGrid)

   call calcderivatives(u, dudt)
end subroutine f_calcderivatives

subroutine f_calcderivativesSquid(u, dudt)
   use FEISTY !, only:
   use globals

   real(dp), intent(in):: u(nGrid)
   real(dp), intent(out):: dudt(nGrid)

   call calcderivativesSquid(u, dudt)
end subroutine f_calcderivativesSquid

subroutine f_getrates(u, dudt, flvl_r, mortpred_r, g_r)
use FEISTY
use globals

         real(dp),intent(in) :: u(nGrid)
         real(dp),intent(inout):: dudt(nGrid)
         real(dp),intent(out) :: flvl_r(nGrid)
         real(dp),intent(out) :: mortpred_r(nGrid)
         real(dp),intent(out) :: g_r(nGrid-nResources)

         call getrates(u, dudt, flvl_r, mortpred_r, g_r)

end subroutine f_getrates

