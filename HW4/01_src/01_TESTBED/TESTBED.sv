/******************************************************************************
* Copyright (C) 2026 Marco
*
* File Name:    TESTBED.sv
* Project:      NTU-DCIC HW4 Symbol Boundary Detector
* Module:       TESTBED
* Author:       Marco <harry2963753@gmail.com>
*
******************************************************************************/
`timescale 1ns/1ps
`include "define.vh"

module TESTBED;

    //=============================================================
    //                      DUT Connections
    //=============================================================
    logic clk, rst_n;
    logic in_valid, in_last;
    logic signed [`W_IN-1:0] in_re, in_im;
    logic [`W_PHI-1:0] eta;
    logic phi_valid;
    logic [`W_IDX-1:0] phi_idx;
    logic [`W_PHI-1:0] phi;
    logic phi_above;
    logic det_valid;
    logic [`W_IDX-1:0] det_idx;
    logic [`W_PHI-1:0] det_phi;

    //=============================================================
    //                         Waveform
    //=============================================================
    initial begin
        $fsdbDumpfile("TESTBED.fsdb");
        $fsdbDumpvars(0, TESTBED, "+mda");
    end

    //=============================================================
    //                   Sim Mode & SDF Annotate
    //=============================================================
    `ifdef GATE
        initial begin
            $display("======================================");
            $display("  [INFO] GATE-LEVEL SIMULATION START  ");
            $display("======================================");
            $sdf_annotate("../02_SYN/Netlist/SBD.sdf", u_dut, , ,"maximum");
        end
    `else
        initial begin
            $display("======================================");
            $display("  [INFO] BEHAVIORAL SIMULATION START  ");
            $display("======================================");
        end
    `endif

    //=============================================================
    //                         Instances
    //=============================================================
    PATTERN u_pattern (
        .clk(clk),
        .rst_n(rst_n),
        .in_valid(in_valid),
        .in_last(in_last),
        .in_re(in_re),
        .in_im(in_im),
        .eta(eta),
        .phi_valid(phi_valid),
        .phi_idx(phi_idx),
        .phi(phi),
        .phi_above(phi_above),
        .det_valid(det_valid),
        .det_idx(det_idx),
        .det_phi(det_phi)
    );

    SBD u_dut (
        .clk(clk),
        .rst_n(rst_n),
        .in_valid(in_valid),
        .in_last(in_last),
        .in_re(in_re),
        .in_im(in_im),
        .eta(eta),
        .phi_valid(phi_valid),
        .phi_idx(phi_idx),
        .phi(phi),
        .phi_above(phi_above),
        .det_valid(det_valid),
        .det_idx(det_idx),
        .det_phi(det_phi)
    );

endmodule
