/** \brief write model run statistics to file
 *
 * \author Magnus Hagdorn
 * \date April 2009
 *
 * \param resname name of the output result file
 * \param cfgname name of the model configuration file
 * \param time    the elapsed wall clock time in seconds
 *
 * open statistics file (create it if it does not exist), gather info
 * from environment and write to file.
 * The file gets locked to avoid parallel access
 */
void gc_writestats(const char *resname, const char *cfgname, double wallTime);
