/*
 * Copyright (c) 2016      Intel, Inc. All rights reserved
 * Copyright (c) 2016      Research Organization for Information Science
 *                         and Technology (RIST). All rights reserved.
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER$
 */

#ifndef PMIX_RENAME_H
#define PMIX_RENAME_H

#define PMIx_Init                      PMIx_Init
#define PMIx_Initialized               PMIx_Initialized
#define PMIx_Finalize                  PMIx_Finalize
#define PMIx_Abort                     PMIx_Abort
#define PMIx_Put                       PMIx_Put
#define PMIx_Commit                    PMIx_Commit
#define PMIx_Fence                     PMIx_Fence
#define PMIx_Fence_nb                  PMIx_Fence_nb
#define PMIx_Get                       PMIx_Get
#define PMIx_Get_nb                    PMIx_Get_nb
#define PMIx_Publish                   PMIx_Publish
#define PMIx_Publish_nb                PMIx_Publish_nb
#define PMIx_Lookup                    PMIx_Lookup
#define PMIx_Lookup_nb                 PMIx_Lookup_nb
#define PMIx_Unpublish                 PMIx_Unpublish
#define PMIx_Unpublish_nb              PMIx_Unpublish_nb
#define PMIx_Spawn                     PMIx_Spawn
#define PMIx_Spawn_nb                  PMIx_Spawn_nb
#define PMIx_Connect                   PMIx_Connect
#define PMIx_Connect_nb                PMIx_Connect_nb
#define PMIx_Disconnect                PMIx_Disconnect
#define PMIx_Disconnect_nb             PMIx_Disconnect_nb
#define PMIx_Resolve_peers             PMIx_Resolve_peers
#define PMIx_Resolve_nodes             PMIx_Resolve_nodes
#define PMIx_Query_info_nb             PMIx_Query_info_nb
#define PMIx_Log_nb                    PMIx_Log_nb

#define PMIx_server_init               PMIx_server_init
#define PMIx_server_finalize           PMIx_server_finalize
#define PMIx_generate_regex            PMIx_generate_regex
#define PMIx_generate_ppn              PMIx_generate_ppn
#define PMIx_server_register_nspace    PMIx_server_register_nspace
#define PMIx_server_deregister_nspace  PMIx_server_deregister_nspace
#define PMIx_server_register_client    PMIx_server_register_client
#define PMIx_server_deregister_client  PMIx_server_deregister_client
#define PMIx_server_setup_fork         PMIx_server_setup_fork
#define PMIx_server_dmodex_request     PMIx_server_dmodex_request

#define PMIx_tool_init                 PMIx_tool_init
#define PMIx_tool_finalize             PMIx_tool_finalize

#define PMIx_Register_event_handler    PMIx_Register_event_handler
#define PMIx_Deregister_event_handler  PMIx_Deregister_event_handler
#define PMIx_Notify_event              PMIx_Notify_event
#define PMIx_Error_string              PMIx_Error_string
#define PMIx_Proc_state_string         PMIx_Proc_state_string
#define PMIx_Persistence_string        PMIx_Persistence_string
#define PMIx_Data_range_string         PMIx_Data_range_string
#define PMIx_Info_directives_string    PMIx_Info_directives_string
#define PMIx_Data_type_string          PMIx_Data_type_string
#define PMIx_Get_version               PMIx_Get_version
#define PMIx_Store_internal            PMIx_Store_internal

#define pmix_value_load                pmix_value_load
#define pmix_value_xfer                pmix_value_xfer
#define pmix_globals                   pmix_globals
#define pmix_output                    pmix_output
#define pmix_output_verbose            pmix_output_verbose

#endif
