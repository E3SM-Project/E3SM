# OMEGA doc/design directory

This directory contains requirement and design documents for the
Ocean Model for E3SM Global Applications (OMEGA) source code.  A
file called OmegaDesignTemplate.rst is provided as a template
for the information required in such a design document.

```mermaid
gantt
    title Omega1 to Omega2 integration plan
    dateFormat  YYYY-MM-DD
    axisFormat  %m/%d
    excludes    weekends

    section Code integration milestones
    INT-000 Omega1                    :milestone, omega1, 2026-04-14, 0d
    INT-010 Surf restoring complete   :milestone, surf_done, 2026-04-15, 0d
    INT-020 Centered pgrad complete   :milestone, pgrad_done, 2026-04-15, 0d
    INT-030 SSH fix integrated        :milestone, ssh_done, 2026-04-24, 0d
    INT-040 Const EOS integrated      :milestone, eos_done, 2026-04-24, 0d
    INT-050 Timestep-EOS integrated  :milestone, tt_done, 2026-05-01, 0d
    INT-055 Timestep-Tend integrated  :milestone, tt_done, 2026-05-14, 0d
    INT-060 Prescribe vel integrated  :milestone, pv_done, 2026-05-21, 0d
    INT-070 FCT integrated            :milestone, fct_done, 2026-05-28, 0d
    INT-080 KPP integrated            :milestone, kpp_done, 2026-05-28, 0d
    INT-090 Submeso eddy integrated   :milestone, submeso_done, 2026-05-28, 0d
    INT-100 Split-explicit integrated :milestone, split_done, 2026-06-04, 0d
    INT-110 Forcing class integrated  :milestone, forcing_done, 2026-06-04, 0d
    INT-120 Surface mass fluxes       :milestone, mass_done, 2026-06-12, 0d
    INT-130 Surface heat fluxes       :milestone, heat_done, 2026-06-12, 0d
    INT-140 Surface salt fluxes       :milestone, salt_done, 2026-06-12, 0d
    INT-150 Partial AM integrated     :milestone, pam_done, 2026-06-23, 0d
    INT-160 Frazil integrated         :milestone, frazil_done, 2026-06-30, 0d
    INT-170 Full AM integrated        :milestone, fam_done, 2026-07-03, 0d
    INT-999 Omega2                    :milestone, omega2, 2026-07-31, 0d

    section Tests of integrated code
    VER-101 Prep RK2 test suite                     :ver101, 2026-04-15, 3d
    VER-102 Prep FB test suite                      :ver102, 2026-04-15, 3d
    VER-103 Run RK2                                 :ver103, after ver101, 4d
    VER-104 Run RK4                                 :ver104, after ver102, 4d

    DEV-401/VER-401 Sphere transport prescribed velocity :devver401, 2026-04-15, 8d

    VER-201 Prep single column w surface restoring  :ver201, after omega1, 4d

    section Omega-Dev1
    DEV-301/VER-301 Implement and test SSH fix      :devver301, after pgrad_done, 7d
    DEV-302/VER-301 Implement and test constant EOS :devver302, after pgrad_done, 7d

    DEV-101/VER-402 Implement and test Timestep-Tend:devver101, after pgrad_done, 9d

    section Omega-Dev2
    DEV-501/VER-501 Implement and test FCT          :devver501, after pv_done, 9d

    DEV-601/VER-601 Implement and test KPP          :devver601, after tt_done, 9d
    section Omega-Dev3 surface fluxes
    DEV-600 Implement Forcing class                 :dev600, after tt_done, 4d
    DEV-602/VER-602 Implement and test Mass fluxes  :devver602, after dev600, 7d
    DEV-603/VER-603 Implement and test Heat fluxes  :devver603, after dev600, 7d
    DEV-604/VER-604 Implement and test Salt fluxes  :devver604, after dev600, 7d

    section Global ocean standalone
    DEV-605/VER-605 Implement and test HO HPGF      :devver605, 2026-05-18, 9d
    VER-510 HR Global wind-forced                   :ver510, after devver605, 6d
    VER-630 HR Global all-forcing                   :ver630, after ver510, 7d

    VER-610 LR IC conversion                        :ver610, 2026-06-01, 5d
    VER-620 LR global simulation                    :ver620, after split_done, 7d

    section Global ocean coupled analysis members
    DEV-701 Time averaging                          :dev701, after split_done, 5d
    DEV-705 Flux accum                              :dev705, after split_done, 5d
    DEV-702 Conservation                            :dev702, after pam_done, 5d
    DEV-703 MOC                                     :dev703, after pam_done, 5d
    DEV-704 Eddy                                    :dev704, after pam_done, 5d

    section Global ocean coupled frazil
    DEV-606/VER-606 Implement and test frazil       :devver606, after tt_done, 9d

    section Coupler Dev1
    DEV-801 Implement MCT momentum                  :dev801, after pam_done, 5d
    DEV-802 Implement full MCT                      :dev802, after dev801, 5d
    DEV-803 Implement Coupler budgets               :dev803, after dev802, 4d

    section Coupled simulation
    VER-700 Coupled wind-forced                     :ver700, after dev801, 7d
    VER-701 HR Coupled all-forced                   :ver701, after dev802, 8d
```

```mermaid
flowchart LR
    classDef VER fill:#99CCFF,stroke-width:2px;
    classDef INT fill:#CCFFCC,stroke-width:2px;
    classDef DEV fill:#FCCCCC,stroke-width:2px;

    subgraph CI ["Code integration"]
        style CI fill:#FFF,stroke:#333,stroke-width:2px
        A@{ shape: circ, label: "Omega1" } ---> B@{ shape: circ, label: "Surf restoring ✔" }
        B ---> C@{ shape: circ, label: "Centered pgrad ✔" }
        C --> E@{ shape: circ, label: "SSH fix" }
        C --> D@{ shape: circ, label: "Const EOS" }
        C --> F@{ shape: circ, label: "Timestep-Tend" }
        F --> G@{ shape: circ, label: "Prescribe vel" }
        G --> H@{ shape: circ, label: "FCT" }
        F --> I@{ shape: circ, label: "KPP" }
        I --> K
        F --> J@{ shape: circ, label: "Submeso eddy" }
        J --> K@{ shape: circ, label: "Split-explicit" }
        K --> Z
        F --> P@{ shape: circ, label: "Forcing class" }
        P --> L@{ shape: circ, label: "Surface mass fluxes" }
        P --> M@{ shape: circ, label: "Surface heat fluxes" }
        P --> N@{ shape: circ, label: "Surface salt fluxes" }
        L --> U@{ shape: circ, label: "Partial AM" }
        M --> U
        N --> U
        U --> W@{ shape: circ, label: "Frazil" }
        W --> V
        U --> V@{ shape: circ, label: "Full AM" }
        V --> Z@{ shape: circ, label: "Omega2" }
        class A,B,C,D,E,F,G,H,I,J,K,L,M,N,P,U,V,W,Z INT
    end

    subgraph TI ["Tests of integrated code"]
        style TI fill:#FFF,stroke:#333,stroke-width:2px
        subgraph Tracer advection STD
            A --> DEV401[Sphere transport prescribed velocity]
            DEV401 --> VER401[Run sphere transport]
            VER401 --> G
            class DEV401 DEV
            class VER401 VER
            end
        subgraph Timestepper
            A --> VER101[Prep RK2 test suite]
            A --> VER102[Prep FB test suite]
            VER101 --> VER103[Run RK2]
            VER102 --> VER104[Run RK4]
            class VER101,VER102,VER103,VER104 VER
            end
        subgraph Surface restoring
            B --> VER201[Prep single column w surface restoring]
            VER201 --> VER202(Run single column)
            class VER201,VER202 VER
            end
    end

    subgraph O1 ["Omega-DEV1"]
        style O1 fill:#FFF,stroke:#333,stroke-width:2px
        subgraph Timestep-Tend
            C --> DEV101[Implement Timestep-Tend]
            DEV101 --> VER402(Run overflow, seamount)
            VER402 --> F
            class DEV101 DEV
            class VER402 VER
            end

        subgraph SSH fix
            DEV302[Implement constant EOS ✔] --> VER301
            C --> DEV301[Implement SSH fix ✔]
            DEV301 --> VER301[Run barotropic channel]
            VER301 --> D
            VER301 --> E
            class DEV301,DEV302 DEV
            class VER301 VER
            end
    end

    subgraph O2 ["Omega-Dev2"]
        style O2 fill:#FFF,stroke:#333,stroke-width:2px
        subgraph Tracer advection FCT
            DEV501[Implement FCT] --> VER501[Run sphere transport w FCT]
            G --> VER501
            VER501 --> H
            class DEV501 DEV
            class VER501 VER
        end
        subgraph KPP
            DEV601[Implement KPP] --> VER601[Run single column w KPP]
            VER601 --> I
            class DEV601 DEV
            class VER601 VER
        end
    end

    subgraph GL ["Global ocean standalone"]
        style GL fill:#FFF,stroke:#333,stroke-width:2px
        subgraph BS ["Bespoke Dev1"]
            VER610[LR IC converstion] --> VER620
            K --> VER620[LR global simulation]
            class VER610,VER620 VER
            end
        subgraph HP ["HPGF"]
            DEV605[Implement HPGF] --> VER605[2-col HPGF test]
            class DEV605 DEV
            class VER605 VER
            end
        subgraph SW ["Global Ocean Mom-only"]
            VER510[HR Global wind-forced]
            VER605 --> VER510
            class VER510 VER
        end
        subgraph PO ["Polaris Dev1"]
            VER630[HR Global all-forcing]
            VER510 --> VER630
            class VER630 VER
            end
    end

    subgraph GC ["Global ocean coupled"]
        style GC fill:#FFF,stroke:#333,stroke-width:2px
        subgraph AM ["Analysis members"]
            K ~~~ DEV701
            DEV701[Time averaging] --> U
            DEV702[Conservation] --> V
            DEV703[MOC] --> V
            DEV704[Eddy] --> V
            DEV705[Flux accum] --> U
            class DEV701,DEV702,DEV703,DEV704,DEV705 DEV
            end
        subgraph CP ["Coupled simulation"]
            P --> VER700[Coupled wind-forced]
            U --> VER701[HR Coupled all-forced]
            V --> VER701
            VER701 ~~~ Z
            class VER700,VER701 VER
        end
        subgraph FR ["Frazil"]
            F --> DEV606[Implement frazil]
            DEV606 --> VER606[Run single column]
            VER606 --> Z
            class DEV606 DEV
            class VER606 VER
        end
        subgraph CPL ["Coupler Dev1"]
            U --> DEV801[Implement MCT momentum]
            DEV801 --> VER700
            DEV801 --> DEV802[Implement full MCT]
            DEV802 --> VER701
            DEV802 --> DEV803[Implement Coupler budgets]
            DEV803 --> VER701
            class DEV801,DEV802,DEV803 DEV
        end
    end

    subgraph O3 ["Omega-Dev3"]
        style O3 fill:#FFF,stroke:#333,stroke-width:2px
        subgraph SF ["Surface fluxes"]
            F --> DEV600[Implement Forcing class]
            DEV600 --> P
            P --> DEV602[Implement Mass fluxes]
            P --> DEV603[Implement Heat fluxes]
            P --> DEV604[Implement Salt fluxes]
            DEV602 --> VER602[Run Single column]
            DEV603 --> VER603[Run Single column]
            DEV604 --> VER604[Run Single column]
            VER602 --> L
            VER603 --> M
            VER604 --> N
            class DEV600,DEV602,DEV603,DEV604 DEV
            class VER602,VER603,VER604 VER
        end
    end
```
