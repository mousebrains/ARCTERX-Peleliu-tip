#! /usr/bin/env bash
# Re-attempt the reference downloads.
#
# Most of these sit behind a subscription, so this will only fetch them from a
# machine with institutional access (campus network or VPN).  Run it again from
# there and the missing PDFs should land next to this script.
#
# It does NOT spoof user agents or route around paywalls.  A 403 here usually
# means the publisher blocked an automated request, not that you lack access --
# open the URL in a browser and the download normally works.
#
#   ./fetch_papers.sh          try everything that is missing
#   ./fetch_papers.sh --list   just print what is missing

set -uo pipefail
cd "$(dirname "$0")" || exit 1

# name|url|citation
ENTRIES=(
"Weiss1991_PhysicaD_enstrophy_transfer.pdf|https://courses.physics.ucsd.edu/2017/Winter/physics216/Weiss-Enstrophy_Transfer.pdf|Weiss 1991, Physica D 48(2-3), 273-294"
"OkuboEbbesmeyer1976_DSR_drogue_vorticity.pdf|https://doi.org/10.1016/0011-7471(76)90875-5|Okubo & Ebbesmeyer 1976, Deep-Sea Res 23(4), 349-352"
"MolinariKirwan1975_JPO_differential_kinematics.pdf|https://doi.org/10.1175/1520-0485(1975)005<0483:CODKPF>2.0.CO;2|Molinari & Kirwan 1975, J Phys Oceanogr 5(3), 483-491"
"Okubo1970_DSR_velocity_singularities.pdf|https://doi.org/10.1016/0011-7471(70)90059-8|Okubo 1970, Deep-Sea Res 17(3), 445-454"
"EfronGong1983_AmStat_bootstrap_jackknife.pdf|https://doi.org/10.1080/00031305.1983.10483087|Efron & Gong 1983, Am Statistician 37(1), 36-48"
# Open access (CC BY); this one should fetch from anywhere, no subscription.
"Poulain2023_FrontMarSci_cyprus_gyre_wavelet_ridge.pdf|https://www.frontiersin.org/articles/10.3389/fmars.2023.1266040/pdf|Poulain et al. 2023, Front Mar Sci 10, 1266040"
# Oceanography 32(4), the FLEAT special issue -- all CC BY 4.0, fetch anywhere.
"Johnston2019_Oceanography_peleliu_wake_eddies_lee_waves.pdf|https://tos.org/oceanography/assets/docs/32-4_johnston2.pdf|Johnston et al. 2019, Oceanography 32(4), 110-125 -- OUR SITE"
"Rudnick2019_Oceanography_vorticity_flow_past_island.pdf|https://tos.org/oceanography/assets/docs/32-4_rudnick.pdf|Rudnick et al. 2019, Oceanography 32(4), 66-73"
"StLaurent2019_Oceanography_palau_wake_turbulence_vorticity.pdf|https://tos.org/oceanography/assets/docs/32-4_st-laurent.pdf|St. Laurent et al. 2019, Oceanography 32(4), 102-109"
"Siegelman2019_Oceanography_palau_near_inertial_surface.pdf|https://tos.org/oceanography/assets/docs/32-4_siegelman.pdf|Siegelman et al. 2019, Oceanography 32(4), 74-83"
"Johnston2019_Oceanography_FLEAT_program_overview.pdf|https://tos.org/oceanography/assets/docs/32-4_johnston1.pdf|Johnston et al. 2019, Oceanography 32(4), 10-21 -- FLEAT overview"
# NOAA Institutional Repository copy.
"Essink2022_JTECH_drifter_cluster_kinematics.pdf|https://repository.library.noaa.gov/view/noaa/65038/noaa_65038_DS1.pdf|Essink et al. 2022, J Atmos Ocean Tech 39(8), 1183-1198"
# Author-posted copy (Univ. of Hawai'i).
"KloosterzielVanHeijst1991_JFM_unstable_barotropic_vortices.pdf|http://www.soest.hawaii.edu/oceanography/faculty/kloosterziel/pdfs/unstable.pdf|Kloosterziel & van Heijst 1991, J Fluid Mech 223, 1-24"
)

# Wanted but publisher-blocked (HTTP 403 to curl); open these in a browser on a
# machine with institutional access.  See "Wanted" in README.md.
#   Zeiden et al. 2022,  J Phys Oceanogr 52(9) 2237-2255   10.1175/JPO-D-21-0252.1
#   Huntley et al. 2022, J Atmos Ocean Tech 39(10) 1499-1523  10.1175/JTECH-D-21-0161.1
#   Ohlmann et al. 2017, Geophys Res Lett 44(1) 330-337    10.1002/2016GL071537
#   MacKinnon et al. 2019, J Geophys Res Oceans 124(7) 4891-4903  10.1029/2019JC014945
#   Zeiden et al. 2019,  J Phys Oceanogr 49(9) 2217-2235   10.1175/JPO-D-18-0233.1
#   Siegelman et al. 2023, J Phys Oceanogr 53(2) 433-455   10.1175/JPO-D-21-0310.1

missing=0
for entry in "${ENTRIES[@]}"; do
    IFS='|' read -r name url cite <<< "$entry"
    if [[ -f "$name" ]]; then
        printf '  have    %s\n' "$name"
        continue
    fi
    if [[ "${1:-}" == "--list" ]]; then
        printf '  MISSING %s\n          %s\n          %s\n' "$name" "$cite" "$url"
        missing=$((missing + 1))
        continue
    fi

    printf '  fetching %s ... ' "$name"
    code=$(curl -sL --max-time 60 -o "$name.part" -w '%{http_code}' "$url" 2>/dev/null)
    if [[ "$code" == "200" ]] && [[ "$(file -b --mime-type "$name.part")" == "application/pdf" ]]; then
        mv "$name.part" "$name"
        printf 'ok\n'
    else
        rm -f "$name.part"
        printf 'not available (HTTP %s)\n           %s\n           %s\n' "$code" "$cite" "$url"
        missing=$((missing + 1))
    fi
done

echo
if (( missing )); then
    echo "$missing reference(s) still missing -- see README.md for the copy-paste list."
else
    echo "All references present."
fi

# Saffman (1992) Vortex Dynamics is a book (ISBN 978-0-521-42058-7); no download.
