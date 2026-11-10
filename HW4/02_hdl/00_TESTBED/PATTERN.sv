/******************************************************************************
* Copyright (C) 2026 Marco
*
* File Name:    PATTERN.sv
* Project:      NTU-DCIC HW4 Symbol Boundary Detector
* Module:       PATTERN
* Author:       Marco <harry2963753@gmail.com>
*
******************************************************************************/
`timescale 1ns/1ps
`include "define.vh"

module PATTERN (
    output logic clk,
    output logic rst_n,
    output logic in_valid,
    output logic in_last,
    output logic signed [`W_IN-1:0] in_re,
    output logic signed [`W_IN-1:0] in_im,
    output logic [`W_PHI-1:0] eta,
    input phi_valid,
    input [`W_IDX-1:0] phi_idx,
    input [`W_PHI-1:0] phi,
    input phi_above,
    input det_valid,
    input [`W_IDX-1:0] det_idx,
    input [`W_PHI-1:0] det_phi
);

    //=============================================================
    //                         Parameters
    //=============================================================
    localparam CLK_PERIOD = 10;
    localparam TIMEOUT = 100000;
    localparam NUM_SMP = 480;
    localparam MAX_DET = 16;
    localparam string VEC_DIR = "../00_TESTBED/src/";

    //=============================================================
    //                       Golden Vectors
    //=============================================================
    logic [2*`W_IN-1:0] in_mem [0:NUM_SMP-1];
    logic [`W_PHI-1:0] phi_mem [0:NUM_SMP-1];
    logic [`W_PHI-1:0] eta_mem [0:0];
    int gold_det_idx [0:MAX_DET-1];
    int gold_det_phi [0:MAX_DET-1];
    int gold_det_num;

    //=============================================================
    //                         Bookkeeping
    //=============================================================
    typedef enum int {PHASE_CLEAN, PHASE_NOISY, DONE} phase_type;
    phase_type veri_phase;
    int cycle, phi_cnt, det_cnt;
    int in_cycle [1:NUM_SMP];
    int out_cycle [1:NUM_SMP];
    int det_cycle [0:MAX_DET-1];
    logic checking;

    //=============================================================
    //                           Clock
    //=============================================================
    initial clk = 1'b0;
    always #(CLK_PERIOD / 2.0) clk = ~clk;

    always @(posedge clk) cycle <= cycle + 1;

    //=============================================================
    //                          Watchdog
    //=============================================================
    initial begin
        #(CLK_PERIOD * TIMEOUT);
        $fatal(1, "[ERROR]: simulation timeout after %0d cycles.", TIMEOUT);
    end

    //=============================================================
    //                         Main Flow
    //=============================================================
    initial begin
        cycle = 0;
        checking = 1'b0;
        veri_phase = PHASE_CLEAN;
        run_case("clean");
        $display("Phase 1 (noiseless) Pass.");
        veri_phase = PHASE_NOISY;
        run_case("noisy");
        $display("Phase 2 (SNR 12 dB) Pass.");
        veri_phase = DONE;
        end_task();
    end

    //=============================================================
    //                           Tasks
    //=============================================================
    task automatic drive_reset();
        in_valid = 1'b0;
        in_last = 1'b0;
        in_re = '0;
        in_im = '0;

        rst_n = 1'b1;
        force clk = 1'b0;
        #20 rst_n = 1'b0;
        #20 rst_n = 1'b1;
        release clk;
        @(negedge clk);
    endtask

    task automatic load_case(input string name);
        int fd, n;
        $readmemh({VEC_DIR, name, "_in.hex"}, in_mem);
        $readmemh({VEC_DIR, name, "_phi.hex"}, phi_mem);
        $readmemh({VEC_DIR, "eta.hex"}, eta_mem);
        fd = $fopen({VEC_DIR, name, "_det.txt"}, "r");
        if(fd == 0) $fatal(1, "[ERROR]: cannot open %s_det.txt", name);
        n = $fscanf(fd, "%d", gold_det_num);
        for(int k = 0; k < gold_det_num; k++) n = $fscanf(fd, "%d %d", gold_det_idx[k], gold_det_phi[k]);
        $fclose(fd);
    endtask

    task automatic run_case(input string name);
        int t0;
        load_case(name);
        eta = eta_mem[0];
        phi_cnt = 0;
        det_cnt = 0;
        drive_reset();
        checking = 1'b1;
        t0 = cycle;
        for(int i = 0; i < NUM_SMP; i++) begin
            in_valid = 1'b1;
            in_last = (i == NUM_SMP - 1);
            in_re = in_mem[i][2*`W_IN-1:`W_IN];
            in_im = in_mem[i][`W_IN-1:0];
            in_cycle[i+1] = cycle;
            @(negedge clk);
        end
        in_valid = 1'b0;
        in_last = 1'b0;
        in_re = '0;
        in_im = '0;
        while(phi_cnt < NUM_SMP || det_cnt < gold_det_num) @(negedge clk);
        repeat(5) @(negedge clk);
        checking = 1'b0;
        if(det_cnt != gold_det_num) $fatal(1, "[ERROR]: %s: %0d detections, golden %0d", name, det_cnt, gold_det_num);
        $display("[%s] eta = %0d, %0d metric samples and %0d detections match the bit-true model", name, eta, phi_cnt, det_cnt);
        $display("[%s] metric latency (input sample -> phi_valid) = %0d cycles", name, out_cycle[1] - in_cycle[1]);
        for(int k = 0; k < det_cnt; k++) begin
            $display("[%s] boundary %0d: m = %0d, Phi_NM = %0d/2^14, det_valid %0d cycles after sample m entered, %0d cycles after its phi_valid", name, k + 1,
                gold_det_idx[k], gold_det_phi[k], det_cycle[k] - in_cycle[gold_det_idx[k]], det_cycle[k] - out_cycle[gold_det_idx[k]]);
        end
    endtask

    task automatic end_task();
        $display("======================================");
        $display("        Congratulations! ALL PASS     ");
        $display("======================================");
        $finish;
    endtask

    //=============================================================
    //                        Output Check
    //=============================================================
    always @(negedge clk) begin
        if(checking) begin
            if(phi_valid) begin
                phi_cnt = phi_cnt + 1;
                out_cycle[phi_idx] = cycle;
                if(phi_idx != phi_cnt) $fatal(1, "[ERROR]: phi_idx = %0d, expected %0d", phi_idx, phi_cnt);
                if(phi !== phi_mem[phi_idx-1]) $fatal(1, "[ERROR]: m = %0d: phi = %0d, golden %0d", phi_idx, phi, phi_mem[phi_idx-1]);
                if(phi_above !== (phi > eta)) $fatal(1, "[ERROR]: m = %0d: phi_above = %0b", phi_idx, phi_above);
            end
            if(det_valid) begin
                if(det_cnt >= gold_det_num) $fatal(1, "[ERROR]: extra detection at m = %0d", det_idx);
                if(det_idx != gold_det_idx[det_cnt] || det_phi != gold_det_phi[det_cnt]) begin
                    $fatal(1, "[ERROR]: detection %0d: (m, phi) = (%0d, %0d), golden (%0d, %0d)", det_cnt + 1, det_idx, det_phi,
                        gold_det_idx[det_cnt], gold_det_phi[det_cnt]);
                end
                det_cycle[det_cnt] = cycle;
                det_cnt = det_cnt + 1;
            end
        end
    end

    //=============================================================
    //                 SystemVerilog Assertion
    //=============================================================
    // Property Name: assert property(condition) <pass event> else <fail event>

    OUT_LOW_IN_RESET: assert property(
        @(posedge clk)
        !rst_n |-> (!phi_valid && !det_valid)
    )
    else $fatal(1, "[ERROR]: output valid asserted during reset.");

    DET_AFTER_PHI: assert property(
        @(posedge clk) disable iff(!rst_n)
        det_valid |-> $past(phi_valid)
    )
    else $fatal(1, "[ERROR]: det_valid must follow a metric sample.");

    DET_SINGLE_CYCLE: assert property(
        @(posedge clk) disable iff(!rst_n)
        det_valid |=> !det_valid
    )
    else $fatal(1, "[ERROR]: det_valid must be a one-cycle pulse.");

endmodule
