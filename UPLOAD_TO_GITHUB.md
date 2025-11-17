# Upload to GitHub Instructions

## Step 1: Create GitHub Repository

1. Go to https://github.com/new
2. Repository name: `clb-detection-toolkit`
3. Description: `Comprehensive bioinformatics pipelines for detecting and characterizing colibactin (CLB) biosynthesis genes in metagenomic samples`
4. **Public** (or Private if you prefer)
5. **DO NOT** initialize with README, .gitignore, or license (we already have these)
6. Click "Create repository"

## Step 2: Push Local Repository to GitHub

After creating the repository, GitHub will show you instructions. Use these commands:

```bash
# Navigate to the repository
cd /tmp/clb-detection-toolkit

# Add GitHub remote (replace YOUR_USERNAME with your GitHub username)
git remote add origin https://github.com/YOUR_USERNAME/clb-detection-toolkit.git

# Rename branch to main (optional but recommended)
git branch -M main

# Push to GitHub
git push -u origin main
```

**Alternative using SSH (if you have SSH keys configured):**

```bash
git remote add origin git@github.com:YOUR_USERNAME/clb-detection-toolkit.git
git branch -M main
git push -u origin main
```

## Step 3: Verify Upload

Visit your repository at:
```
https://github.com/YOUR_USERNAME/clb-detection-toolkit
```

You should see:
- ✓ README.md displayed on the homepage
- ✓ All pipeline directories (pipeline1_bbmap, pipeline2_bowtie2, pipeline3_assembly)
- ✓ Utils and example configs
- ✓ MIT License

## Optional: Add Topics/Tags

On your GitHub repository page, click "Add topics" and add:
- `bioinformatics`
- `metagenomics`
- `colibactin`
- `microbiome`
- `genomics`
- `pipeline`

This helps others discover your toolkit!

## Troubleshooting

**If you get authentication errors:**

1. **Using HTTPS**: You may need to create a personal access token
   - Go to Settings → Developer settings → Personal access tokens → Generate new token
   - Select scopes: `repo` (full control of private repositories)
   - Use token as password when pushing

2. **Using SSH**: Make sure your SSH keys are set up
   - Check: `ssh -T git@github.com`
   - Setup guide: https://docs.github.com/en/authentication/connecting-to-github-with-ssh

**If you need to change the remote URL:**

```bash
# Check current remote
git remote -v

# Change to HTTPS
git remote set-url origin https://github.com/YOUR_USERNAME/clb-detection-toolkit.git

# Or change to SSH
git remote set-url origin git@github.com:YOUR_USERNAME/clb-detection-toolkit.git
```
