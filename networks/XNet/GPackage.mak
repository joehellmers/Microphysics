XNET_sources += xnet_abundances.F90 \
                xnet_common.F90 \
                xnet_conditions.F90 \
                xnet_constants.F90 \
                xnet_control.F90 \
                xnet_data.F90 \
                xnet_data_distribute_mpi.F90 \
                xnet_eos.F90 \
                xnet_fd.F90 \
                xnet_ffn.F90 \
                xnet_finalize.F90 \
                xnet_flux.F90 \
                xnet_full_net.F90 \
                xnet_init.F90 \
                xnet_interface.F90 \
                xnet_interface_mpi.F90 \
                xnet_jacobian_bcast_dense.F90 \
                xnet_jacobian_dense.F90 \
                xnet_match.F90 \
                xnet_net_preprocess.F90 \
                xnet_output.F90 \
                xnet_screening.F90 \
                xnet_solve_be.F90 \
                xnet_terminate.F90 \
                xnet_types.F90 \
                xnet_util.F90

F90sources += actual_network.F90
ifneq ($(USE_REACT), FALSE)
  F90sources += actual_burner.F90
  F90sources += actual_rhs.F90
endif
F90sources += $(XNET_sources)

# actual_network.F90 is created at build time for this network
actual_network.F90: $(MICROPHYSICS_HOME)/networks/XNet/network.template
	@echo " "
	@echo "---------------------------------------------------------------------------"
	@echo "${bold}WRITING actual_network.F90${normal}"
	$(MICROPHYSICS_HOME)/networks/XNet/write_network.py \
            -t $(MICROPHYSICS_HOME)/networks/XNet/network.template \
            -s $(MICROPHYSICS_HOME)/networks/XNet/Networks/$(XNET_DATA) \
            -o actual_network.F90
	@echo "---------------------------------------------------------------------------"
	@echo " "


# remove actual_network.F90 for 'make clean' and therefore 'make realclean'
clean::
	$(RM) actual_network.F90