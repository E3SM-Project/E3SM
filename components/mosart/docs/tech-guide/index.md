# MOSART Technical Guide

This guide provides scientific and technical details about MOSART.

## Physics

MOSART is a one-dimension river transport model that is designed for river routing at local, regional, and global scales ([Li et al., 2013](https://doi.org/10.1175/JHM-D-12-015.1)). Its primary function is to supply freshwater inputs to ocean models within coupled Earth System Models.

MOSART divides each spatial unit, such as a latitude/longitude grid or a sub-basin, into three hydrologic categories: hillslopes, tributaries, and a main channel (see figure below). The hillslopes receive runoff and send into tributaries, which then converge into a single main channel. This main channel connects adjacent upstream and downstream units through the river network. MOSART simplifies the multiple tributaries within a spatial unit into a single hypothetical sub-network channel, which has a transport capacity equivalent to all combined tributaries.

![alt text](../figures/mosart_concept.png)

- Hillslope Routing: Within each spatial unit, surface runoff is directed as overland flow to the sub-network channel, while subsurface runoff enters the sub-network channel directly.

- Sub-network Channel Routing: This channel aggregates water from the hillslopes, routes it through the channel system, and discharges it into the main channel.

- Main Channel Routing: The main channel collects water from the sub-network channel and any inflow from upstream spatial units, eventually discharging the accumulated water downstream to the next spatial unit or directly to the ocean.

## Parameters

### Main parameters required in the MOSART parameter file

| Variable Name     | Description [unit]                                                             |
|-------------------|--------------------------------------------------------------------------------|
| fdir              | Flow direction [unitless]                                                      |
| lat               | Latitude at cell center [degree]                                               |
| lon               | Longitude at cell center [degree]                                              |
| frac              | fraction of the unit draining to the outlet [0-1]                              |
| rslp              | main channel slope [unitless]                                                  |
| rlen              | main channel length [m]                                                        |
| tslp              | mean tributary channel slope averaged through the unit [unitless]              |
| area              | local drainage area [m^2]                                                      |
| areaTotal         | total upstream drainage area, local unit included; multi flow direction [m^2]  |
| areaTotal2        | total upstream drainage area, local unit included; single flow direction [m^2] |
| rdep              | bankfull depth of main channel [m]                                             |
| rwid              | bankfull width of main channel [m]                                             |
| rwid0             | floodplain width linked to main channel [m]                                    |
| gxr               | drainage density [m^-1]                                                        |
| hslp              | topographic slope [unitless]                                                   |
| twid              | bankfull width of local tributaries [m]                                        |
| nr                | Manning''s roughness coefficient for main channel flow [unitless]              |
| nt                | Manning''s roughness coefficient for tributary channel flow [unitless]         |
| nh                | Manning''s roughness coefficient for overland flow [unitless]                  |

### Parameters required by additional MOSART features

Coming soon.
