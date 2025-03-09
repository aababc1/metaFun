#!/bin/bash

# Set variables
PACKAGE_NAME="metafun_test"
CONDA_ENV_PATH="/data1/leehg/miniconda3/envs/build_env"
BUILD_PATH="${CONDA_ENV_PATH}/conda-bld"

echo "============================================="
echo "Building and uploading metafun_test package"
echo "============================================="

# Step 1: Clean up old builds
echo "Removing old builds..."
rm -rf "${BUILD_PATH}"
if [ $? -eq 0 ]; then
  echo "✓ Old builds removed successfully"
else
  echo "! Warning: Could not remove old builds, continuing anyway"
fi

/data1/leehg/miniconda3/envs/build_env/bin/conda build purge


# Step 2: Build the package
echo "Building package..."
${CONDA_ENV_PATH}/bin/conda build --no-test ./
if [ $? -ne 0 ]; then
  echo "✗ Build failed! Exiting."
  exit 1
fi
echo "✓ Package built successfully"

# Step 3: Find the new package
PACKAGE_PATH=$(find ${BUILD_PATH} -name "${PACKAGE_NAME}*.tar.bz2" | sort -r | head -n 1)
if [ -z "$PACKAGE_PATH" ]; then
  echo "✗ Could not find built package! Exiting."
  exit 1
fi
echo "✓ Found package: $PACKAGE_PATH"

# Step 4: Upload to Anaconda Cloud
echo "Uploading to Anaconda Cloud..."
anaconda upload --force $PACKAGE_PATH
if [ $? -eq 0 ]; then
  echo "✓ Upload successful!"
else
  echo "✗ Upload failed. You may need to run 'anaconda login' first."
  exit 1
fi

echo "============================================="
echo "Build and upload process complete!"
echo "============================================="
