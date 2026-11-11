#!/usr/bin/env bash
# Download and unpack the IonStats example dataset (~200 MB) into ./example/.
set -euo pipefail

# figshare's own API redirects file/download/<id> to this host, so it is the
# canonical form. The figshare.com/ndownloader/... path answers a bot check with
# an empty 202 instead of the file, which used to leave a 0-byte archive here.
#
# Do NOT add --user-agent "Mozilla/5.0". This host redirects to presigned
# Cloudflare R2 storage, and Cloudflare 403s that string: measured on the same
# URL, curl's own user agent returns 206/200 while "Mozilla/5.0" returns 403.
URL="https://ndownloader.figshare.com/files/59204378"
ARCHIVE="example_data.tar.gz"

fetch() {
    if command -v curl > /dev/null 2>&1; then
        curl -fSL --retry 3 --retry-delay 5 -o "$1" "$2"
    elif command -v wget > /dev/null 2>&1; then
        wget --tries=3 -O "$1" "$2"
    else
        echo "ERROR: neither curl nor wget is installed." >&2
        exit 1
    fi
}

echo "Downloading the example dataset (~200 MB) from figshare ..."
fetch "$ARCHIVE" "$URL"

# An error or challenge page can arrive with a success status, which leaves a
# small HTML file named like an archive. Check before unpacking, and do not
# leave a stale file behind for the next run to trip over.
if [ ! -s "$ARCHIVE" ]; then
    echo "ERROR: $ARCHIVE is empty -- the server did not return the file." >&2
    echo "Open $URL in a browser to check that the item is public." >&2
    rm -f "$ARCHIVE"
    exit 1
fi
if [ "$(head -c 2 "$ARCHIVE" | od -An -tx1 | tr -d ' \n')" != "1f8b" ]; then
    echo "ERROR: $ARCHIVE is not a gzip archive. It starts with:" >&2
    head -c 200 "$ARCHIVE" >&2
    echo >&2
    rm -f "$ARCHIVE"
    exit 1
fi

tar -xzf "$ARCHIVE"
rm -f "$ARCHIVE"
echo "Done. The example data is in ./example/"
