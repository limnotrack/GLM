/******************************************************************************
 * glm_lib.c - Library API for GLM (callable from Python, etc.)
 *
 * Provides glm_run_with_nml() so GLM can be invoked as a library without
 * the main() entry point. Used by pglm Python package.
 ******************************************************************************/
#include <string.h>

extern char glm_nml_file[];
extern void run_model(void);

/* Stub for glm_output.c (defined in glm_main.c for executable build) */
char *all_plots_name = NULL;

#ifdef _WIN32
__declspec(dllexport)
#endif
void glm_run_with_nml(const char *nml_path)
{
    if (nml_path && nml_path[0]) {
        strncpy(glm_nml_file, nml_path, 255);
        glm_nml_file[255] = '\0';
    }
    run_model();
}
