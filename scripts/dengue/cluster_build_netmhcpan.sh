#!/usr/bin/env bash
# Build NetMHCpan I + II SIFs on cluster.
#
# This is defensive: tries multiple paths because Peter's existing licensed
# tarballs may be in any of several known locations on the cluster.
#
# Strategy (in order):
#   1. Already a SIF? skip.
#   2. Already a docker image? -> apptainer build from docker-daemon://
#   3. Tarball in container-build-items/? -> docker build then SIF convert
#   4. Tarball in /data/license_archives or /data/mccaffrey/...? -> stage + build
#   5. Otherwise SKIP and log; T-cell phase will be marked unavailable.
#
# Usage on cluster:
#   bash scripts/dengue/cluster_build_netmhcpan.sh

set -uo pipefail

REPO=/data/james/ica-dengue/immunoinformatics_platform-dengue
LOG=/data/james/ica-dengue/logs/netmhcpan_build.log
SIF_PATH=/data/james/ica-dengue/sif_images
mkdir -p "$(dirname "$LOG")" "$SIF_PATH"

log() { printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "$LOG"; }

cd "$REPO"
log "netmhcpan build starting"

# Common search paths for Peter's licensed tarballs
search_paths=(
    "$REPO/container-build-items"
    /data/license_archives
    /data/mccaffrey/license_archives
    /data/mccaffrey/container-build-items
    /data/peter/container-build-items
    /home/pathinformatics/container-build-items
    /data/tools/netmhcpan
)

build_one() {
    local class=$1            # i or ii
    local sif_name=$2         # netmhcpan_i.sif
    local docker_tag=$3       # netmhcpan-i (matches build-net-mhc-pan-i.sh from Peter)
    local container_dir=$4    # containers/netmhcpan_i_4.1
    local tarball_glob=$5     # netMHCpan-4.1*.tar.gz

    log ""
    log "----- NetMHCpan $class -----"

    # Check 1: SIF already exists
    if [[ -f "$SIF_PATH/$sif_name" ]]; then
        log "  $sif_name already exists, skipping build"
        return 0
    fi

    # Check 2: docker image already exists (Peter's prior build)
    if command -v docker >/dev/null 2>&1; then
        if docker images --format "{{.Repository}}" | grep -q "^${docker_tag}$"; then
            log "  found docker image $docker_tag; converting to SIF"
            if apptainer build "$SIF_PATH/$sif_name" "docker-daemon://${docker_tag}:latest" >> "$LOG" 2>&1; then
                log "  $sif_name built from docker daemon"
                return 0
            else
                log "  apptainer build from docker failed; will try other paths"
            fi
        fi
    fi

    # Check 3 + 4: tarball in any search path -> stage to container-build-items, then docker build
    log "  searching for $tarball_glob in known locations"
    found_tarball=""
    for sp in "${search_paths[@]}"; do
        if [[ -d "$sp" ]]; then
            match=$(find "$sp" -maxdepth 2 -name "$tarball_glob" 2>/dev/null | head -1)
            if [[ -n "$match" ]]; then
                found_tarball="$match"
                log "    found: $found_tarball"
                break
            fi
        fi
    done

    if [[ -z "$found_tarball" ]]; then
        log "  no tarball found in: ${search_paths[*]}"
        log "  SKIPPING NetMHCpan $class build; T-cell pipeline phase will be unavailable"
        return 1
    fi

    # Stage tarball and build
    mkdir -p "$REPO/container-build-items"
    cp "$found_tarball" "$REPO/container-build-items/" 2>>"$LOG" || true

    log "  building docker image $docker_tag from $container_dir"
    if docker build -t "$docker_tag:latest" \
            -f "$REPO/$container_dir/Dockerfile" \
            "$REPO" >> "$LOG" 2>&1; then
        log "  docker build OK, converting to SIF"
        if apptainer build "$SIF_PATH/$sif_name" "docker-daemon://${docker_tag}:latest" >> "$LOG" 2>&1; then
            log "  $sif_name built"
            return 0
        else
            log "  apptainer build FAILED"
            return 2
        fi
    else
        log "  docker build FAILED; check $LOG"
        return 3
    fi
}

build_one "I"  "netmhcpan_i.sif"  "netmhcpan-i"  "containers/netmhcpan_i_4.1"  "netMHCpan-4.1*.tar.gz"
build_one "II" "netmhcpan_ii.sif" "netmhcpan-ii" "containers/netmhcpan_ii_4.1" "netMHCIIpan-4.*.tar.gz"

log ""
log "netmhcpan build script complete"
log "SIFs in $SIF_PATH:"
ls -lah "$SIF_PATH"/netmhcpan*.sif 2>/dev/null | tee -a "$LOG"
