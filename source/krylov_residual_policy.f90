! Internal opt-in additive residual policy; not a stable library API.
module krylov_residual_policy
  use iso_fortran_env, only: int64
  implicit none
  private
  public :: residual_policy,valid_residual_policy,requested_residual_limit
  type :: residual_policy
    ! Zero fields are storage defaults, NOT a usable default tolerance policy.
    real(8) :: atol=0d0,rtol=0d0,operator_scale=0d0,origin=0d0
    logical :: scale_provided=.false.
  end type
contains
  logical function finite_policy_value(x) result(ok)
    real(8), intent(in) :: x
    ok=iand(transfer(x,0_int64),int(z'7ff0000000000000',int64))/=int(z'7ff0000000000000',int64)
  end function

  logical function valid_residual_policy(p,reason) result(ok)
    type(residual_policy), intent(in) :: p
    character(*), intent(out) :: reason
    real(8) :: scale,threshold
    ok=.false.;reason='nonfinite residual policy parameter'
    if(.not.finite_policy_value(p%atol) .or. .not.finite_policy_value(p%rtol) .or. &
       .not.finite_policy_value(p%operator_scale) .or. .not.finite_policy_value(p%origin))return
    reason='negative residual tolerance or operator scale'
    if(p%atol<0d0 .or. p%rtol<0d0 .or. p%operator_scale<0d0)return
    reason='at least one positive residual tolerance is required'
    if(p%atol==0d0 .and. p%rtol==0d0)return
    reason='relative residual tolerance requires an explicit H-origin*I scale'
    if(p%rtol>0d0 .and. .not.p%scale_provided)return
    reason='residual threshold overflows for the supplied scale'
    if(.not.requested_residual_limit(p,p%origin,scale,threshold))return
    reason='valid opt-in residual policy';ok=.true.
  end function

  logical function requested_residual_limit(p,value,scale,threshold) result(ok)
    type(residual_policy), intent(in) :: p
    real(8), intent(in) :: value
    real(8), intent(out) :: scale,threshold
    real(8) :: term
    ok=.false.;scale=0d0;threshold=0d0
    ! Parameters are validated once before iteration. Avoid even forming an
    ! unused energy difference for an absolute-only request.
    if(.not.finite_policy_value(value))return
    if(p%rtol==0d0)then
      threshold=p%atol;ok=.true.;return
    end if
    if(value>0d0 .and. p%origin<0d0)then
      if(value>huge(1d0)+p%origin)return
    else if(value<0d0 .and. p%origin>0d0)then
      if(value< -huge(1d0)+p%origin)return
    end if
    scale=max(p%operator_scale,abs(value-p%origin))
    if(.not.finite_policy_value(scale))return
    if(scale>1d0)then
      if(p%rtol>huge(1d0)/scale)return
    end if
    term=p%rtol*scale
    if(term>huge(1d0)-p%atol)return
    threshold=p%atol+term
    ok=finite_policy_value(threshold)
    ! No epsilon floor, automatic relaxation, or observed-scale substitution.
  end function
end module
