# Relaxed trajectories for the runs that yielded accepted binders

These are the post-FastRelax hallucination outputs that ProteinMPNN redesigned
into the accepted designs in `../bc1/*/accepted/`. They are the parent structures
of every binder in this campaign, kept because the 5 GB `designs/` tree did not
survive the cluster shutdown.

- `sweep_af3_guess/` — MD1 block with `predict_initial_guess`; source of the 9 MD1 binders
- `guess_clamp_md1md2_state3/` — MD1-MD2 clamp, TICA state3; source of the 5 clamp binders

The remaining runs' trajectories (md1_block_af3, md1_block_sim, md2/md3 blocks,
clamp state1) produced no accepted designs and were not preserved.
