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
)

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
