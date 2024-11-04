#ifndef SCREAM_SHR_INTERFACE_HPP
#define SCREAM_SHR_INTERFACE_HPP

extern "C"
{

#ifdef SCREAM_CIME_BUILD
int shr_get_iosysid_c2f(int atm_id);
int shr_get_iotype_c2f(int atm_id);
int shr_get_rearranger_c2f(int atm_id);
int shr_get_ioformat_c2f(int atm_id);
#endif

} // extern "C"

#endif // SCREAM_SHR_INTERFACE_HPP
